import os
import sys
import haiku as hk
import jax
import jax.numpy as jnp
import pandas as pd
import numpy as np
import h5py
from nucleotide_transformer.pretrained import get_pretrained_segment_nt_model

# Get the folder path from the user input
folder_path = sys.argv[1]

# Check if batch size and batch number are provided, otherwise set defaults
if len(sys.argv) > 2:
    batch_size = int(sys.argv[2])
else:
    batch_size = len([name for name in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, name))])
    
if len(sys.argv) > 3:
    batch_number = int(sys.argv[3])
else:
    batch_number = 1

batch_start = (batch_size * batch_number) - batch_size
batch_end = batch_size * batch_number

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

# The number of DNA tokens (excluding the CLS token prepended) needs to be divisible by 2 to the power of the number of downsampling blocks, i.e., 4.
max_num_nucleotides = 5000
assert max_num_nucleotides % 4 == 0, (
    "The number of DNA tokens (excluding the CLS token prepended) needs to be divisible by 2 to the power of the number of downsampling blocks, i.e., 4."
)

# Adjust the rescaling factor if necessary
if max_num_nucleotides + 1 > 5001:
    inference_rescaling_factor = (max_num_nucleotides + 1) / 2048
else:
    inference_rescaling_factor = None

# Iterate over each file in the specified folder and process within batch range
for file_name in os.listdir(folder_path)[batch_start:batch_end]:
    # Construct the full file path
    file_path = os.path.join(folder_path, file_name)
    
    # Check if it's a file (and not a directory)
    if os.path.isfile(file_path):
        # Assign the file path to the input_file variable
        input_file = file_path
        
        # Now you can work with input_file
        print(f"Processing file: {input_file}")
        
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
        
        inputfile = pd.read_table(input_file, header=None)
        probabilities = []
        
        # Get data and tokenize it
        sequences = list(inputfile[4])
        
        # Replace sequences with too many 'N's with 'N' * 5000
        for i, seq in enumerate(sequences):
            if seq.count('N') > 50:
                sequences[i] = 'N' * 5000
        
        tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences)]
        tokens_str = [b[0] for b in tokenizer.batch_tokenize(sequences)]
        tokens = jnp.stack([jnp.asarray(tokens_ids, dtype=jnp.int32)] * num_devices, axis=0)
        
        # Infer on the sequence
        outs = apply_fn(parameters, keys, tokens)
        logits = outs["logits"]
        
        # Transform logits into probabilities
        probabilities_unit = jnp.asarray(jax.nn.softmax(logits, axis=-1))[...,-1]
        probabilities.append(probabilities_unit)
        
        output_file = input_file + "SB_probabilities" + ".h5"
        
        # Save the 'probabilities' object to an HDF5 file
        with h5py.File(output_file, 'w') as f:
            for i, sublist in enumerate(probabilities):
                grp = f.create_group(f'group_{i}')
                for j, array in enumerate(sublist):
                    grp.create_dataset(f'array_{j}', data=array)
