import os
import sys
import haiku as hk
import jax
import jax.numpy as jnp
import pandas as pd
import numpy as np
import h5py
from nucleotide_transformer.pretrained import get_pretrained_segment_nt_model
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = ".99"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"

# Define the file path directly
file_path = sys.argv[1]

# Initialize CPU as default JAX device, fall back to CPU if GPU is not available
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

# Set up the model with required parameters
max_num_nucleotides = 5000
assert max_num_nucleotides % 4 == 0, (
    "The number of DNA tokens (excluding the CLS token prepended) needs to be divisible by 2 to the power of the number of downsampling blocks, i.e., 4."
)

# Adjust the rescaling factor if necessary
inference_rescaling_factor = (max_num_nucleotides + 1) / 2048 if max_num_nucleotides + 1 > 5001 else None
parameters, forward_fn, tokenizer, config = get_pretrained_segment_nt_model(
    model_name="segment_nt",
    embeddings_layers_to_save=(29,),
    rescaling_factor=inference_rescaling_factor,
    attention_maps_to_save=((1, 4), (7, 10)),
    max_positions=max_num_nucleotides + 1,
    verbose=False
)

forward_fn = hk.transform(forward_fn)
apply_fn = jax.pmap(forward_fn.apply, axis_name='batch', devices=devices, donate_argnums=(0,))
random_key = jax.random.PRNGKey(seed=0)
keys = jax.random.split(random_key, num_devices)
parameters = jax.device_put_replicated(parameters, devices=devices)

# Open the file for reading and the HDF5 file for writing
output_file = sys.argv[2]
with open(file_path, 'r') as file, h5py.File(output_file, 'w') as f:
    for line_idx, line in enumerate(file):
        print(line_idx)
        # Extract the sequence (assumed to be the 5th column)
        fields = line.strip().split('\t')
        sequence = fields[4] if len(fields) > 4 else 'N' * 5000  # Handle missing data
        
        # Replace sequences with too many 'N's with 'N' * 5000
        if sequence.count('N') > 50:
            sequence = 'N' * 5000
        
        tokens_ids = [b[1] for b in tokenizer.batch_tokenize([sequence])]
        tokens_str = [b[0] for b in tokenizer.batch_tokenize([sequence])]
        tokens = jnp.stack([jnp.asarray(tokens_ids, dtype=jnp.int32)] * num_devices, axis=0)
        
        # Infer on the sequence
        outs = apply_fn(parameters, keys, tokens)
        logits = outs["logits"]
        
        # Transform logits into probabilities
        probabilities_unit = jnp.asarray(jax.nn.softmax(logits, axis=-1))[...,-1]
        
        # Append `probabilities_unit` to the HDF5 file directly
        grp = f.create_group(f'group_{line_idx}')
        for j, array in enumerate(probabilities_unit):
            grp.create_dataset(f'array_{j}', data=array)

        logits.block_until_ready()
        del logits, outs  # Free GPU memory
