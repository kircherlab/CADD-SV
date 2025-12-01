import os
import sys
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np
import h5py
from tqdm import tqdm
from nucleotide_transformer.pretrained import get_pretrained_segment_nt_model

# ----------------------------
# Environment setup
# ----------------------------

os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = ".99"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"

# ----------------------------
# Parse input args
# ----------------------------

file_path = sys.argv[1]
output_file = sys.argv[2]

# ----------------------------
# Device setup
# ----------------------------

try:
    backend = "gpu"
    devices = jax.devices(backend)
    if not devices:
        raise RuntimeError("No GPU devices found.")
except RuntimeError:
    backend = "cpu"
    devices = jax.devices(backend)

num_devices = len(devices)
print(f"Devices found: {devices}")

# ----------------------------
# Model setup
# ----------------------------

max_num_nucleotides = 5000
assert max_num_nucleotides % 4 == 0, "DNA tokens need to be divisible by 4."

inference_rescaling_factor = (
    (max_num_nucleotides + 1) / 2048 if max_num_nucleotides + 1 > 5001 else None
)

parameters, forward_fn, tokenizer, config = get_pretrained_segment_nt_model(
    model_name="segment_nt",
    embeddings_layers_to_save=(29,),
    rescaling_factor=inference_rescaling_factor,
    attention_maps_to_save=((1, 4), (7, 10)),
    max_positions=max_num_nucleotides + 1
)

forward_fn = hk.transform(forward_fn)
apply_fn = jax.pmap(forward_fn.apply, axis_name='batch', devices=devices, donate_argnums=(0,))

random_key = jax.random.PRNGKey(seed=0)
keys = jax.random.split(random_key, num_devices)
parameters = jax.device_put_replicated(parameters, devices=devices)

# ----------------------------
# Processing setup
# ----------------------------

samples_per_device = 8  # ✅ Adjust this for faster processing (try 8, 16, etc.)
batch_size = num_devices * samples_per_device

# ----------------------------
# Count total lines for tqdm
# ----------------------------

with open(file_path) as f:
    total_lines = sum(1 for _ in f)

# ----------------------------
# Main processing loop
# ----------------------------

sequences = []
line_indices = []

with open(file_path, 'r') as infile, h5py.File(output_file, 'w') as f:
    for line_idx, line in enumerate(tqdm(infile, total=total_lines, desc="Processing sequences")):
        fields = line.strip().split('\t')
        sequence = fields[4] if len(fields) > 4 else 'N' * 5000

        # Sanitize sequence
        if sequence.count('N') > 50:
            sequence = 'N' * 5000

        sequences.append(sequence)
        line_indices.append(line_idx)

        # Process batch when full
        if len(sequences) == batch_size:
            tokenized = tokenizer.batch_tokenize(sequences)
            tokens_ids = [b[1] for b in tokenized]
            tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)

            # Reshape for multi-device batch: [num_devices, samples_per_device, seq_length]
            tokens = tokens.reshape((num_devices, samples_per_device, -1))

            outs = apply_fn(parameters, keys, tokens)
            probabilities_unit = jnp.asarray(jax.nn.softmax(outs["logits"], axis=-1))[..., -1]
            probabilities_unit = np.array(probabilities_unit)

            # Write results to HDF5
            idx = 0
            for device_batch in probabilities_unit:
                for prob in device_batch:
                    grp = f.create_group(f'group_{line_indices[idx]}')
                    for j, array in enumerate(prob):
                        grp.create_dataset(f'array_{j}', data=array)
                    idx += 1

            # Cleanup
            logits = outs["logits"]
            logits.block_until_ready()
            del outs, logits, probabilities_unit
            sequences = []
            line_indices = []

    # Process any remaining sequences (if total % batch_size != 0)
    if sequences:
        pad_size = batch_size - len(sequences)
        sequences += ['N' * 5000] * pad_size
        line_indices += [-1] * pad_size  # Dummy index for padded entries

        tokenized = tokenizer.batch_tokenize(sequences)
        tokens_ids = [b[1] for b in tokenized]
        tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)
        tokens = tokens.reshape((num_devices, samples_per_device, -1))

        outs = apply_fn(parameters, keys, tokens)
        probabilities_unit = jnp.asarray(jax.nn.softmax(outs["logits"], axis=-1))[..., -1]
        probabilities_unit = np.array(probabilities_unit)

        # Write results
        idx = 0
        for device_batch in probabilities_unit:
            for prob in device_batch:
                if line_indices[idx] == -1:
                    idx += 1
                    continue  # Skip padded entries
                grp = f.create_group(f'group_{line_indices[idx]}')
                for j, array in enumerate(prob):
                    grp.create_dataset(f'array_{j}', data=array)
                idx += 1

        logits = outs["logits"]
        logits.block_until_ready()
        del outs, logits, probabilities_unit

