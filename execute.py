import os
import shutil
import json

# Configuration settings
system = "CrNi"
undoped = False

# Retrieve size from input file and calculate total atoms
with open("input_file", 'r') as file:
    input_data = {line.split(':')[0].strip(): line.split(':')[1].strip() for line in file if ':' in line}

GB = input_data["grain_boundary"].lower() == 'true'

# Parse additives properly from the input file
if input_data["additives"].lower() == "none":
    additives = []
    undoped = True
else:
    # Convert string format "[(Cr, B, 9), (Ti, H, 9)]" to a list of tuples
    cleaned_value = input_data["additives"].strip()[1:-1]  # Remove brackets
    additives = []
    for additive in cleaned_value.split('), ('):
        components = [v.strip() for v in additive.strip('()').split(',')]
        if len(components) == 3:
            additives.append((components[0], components[1], int(components[2])))
        else:
            raise ValueError(f"Unexpected format for additive: {additive}")

size = int(input_data['size'])
crystal_shape = input_data['crystal_shape']
cuc = {'fcc': 4, 'bcc': 2}
total_atoms = int(size**3 * cuc[crystal_shape])

# Check for POSCAR-gb if GB is True
if GB and not os.path.exists("POSCAR-gb"):
    raise FileNotFoundError("POSCAR-gb file is required but not found. Please ensure the file exists in the working directory.")

dopant_fractions = [0.01, 0.02, 0.04]
if undoped or GB:
    dopant_fractions = [0.01, 0.02, 0.03, 0.04, 0.05]

# Input file path
input_file = "input_file"

# Function to update the input file
def update_input_file(file_path, updates):
    """
    Update the input file with new configuration settings.

    Parameters:
        file_path (str): Path to the input file
        updates (dict): Dictionary of updates to apply

    Returns:
        None
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for key, value in updates.items():
        for i, line in enumerate(lines):
            if line.startswith(key):
                if isinstance(value, dict):
                    value_str = ', '.join([f"{k}={v}" for k, v in value.items()])
                    lines[i] = f"{key}: {value_str}\n"
                elif isinstance(value, list):  # Handle additives list
                    value_str = "[" + ", ".join(f"({x[0]}, {x[1]}, {x[2]})" for x in value) + "]"
                    lines[i] = f"{key}: {value_str}\n"
                elif isinstance(value, tuple):
                    value_str = ', '.join(map(str, value))
                    lines[i] = f"{key}: {value_str}\n"
                else:
                    lines[i] = f"{key}: {value}\n"
                break

    with open(file_path, 'w') as file:
        file.writelines(lines)

# Loop through dopant fractions
for fraction in dopant_fractions:
    num_dopants = round(fraction * total_atoms)

    if GB and not undoped:
        num_dopants = 50

    # Update the additives list with the correct dopant amount
    updated_additives = []
    for metal_host, dopant, _ in additives:
        updated_additives.append((metal_host, dopant, num_dopants))

    # Update the input file
    updates = {
        "composition": {k.strip(): float(v.strip()) for k, v in [item.split('=') for item in input_data['composition'].split(',')]},
        "grain_boundary": GB,
        "additives": updated_additives if not undoped else "None"
    }

    update_input_file(input_file, updates)

    # Execute the simulation
    os.system("python3 hybrid_mcmd.py")

    # Organize results
    dir_name = f"{system}+{'+'.join([dopant for _, dopant, _ in updated_additives])}_{int(fraction * 100)}"
    os.makedirs(dir_name, exist_ok=True)
    for item in ["POSCAR-1", "data", "structures"]:
        if os.path.exists(item):
            if os.path.isdir(item):
                shutil.copytree(item, os.path.join(dir_name, item))
            else:
                shutil.copy(item, os.path.join(dir_name, item))

    print(f"Simulation for {', '.join([dopant for _, dopant, _ in updated_additives])} fraction {fraction*100}% completed.")
