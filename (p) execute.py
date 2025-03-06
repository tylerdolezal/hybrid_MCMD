import os
import shutil
from multiprocessing import Pool, cpu_count
from ase.io import read
from ase.build import bulk

# Read input file
with open("input_file", 'r') as file:
    input_data = {line.split(':')[0].strip(): line.split(':')[1].strip() for line in file if ':' in line}

# Check if batch execution is enabled
batch_fractions = input_data.get("batch_mode", "None").strip()

if batch_fractions.lower() == "none":
    batch_fractions = None
else:
    # Convert string list to actual list of floats
    batch_fractions = [float(x.strip()) for x in batch_fractions.strip("[]").split(",")]

# Extract other parameters
GB = input_data.get("grain_boundary", "False").lower() == "true"
size = int(input_data.get("size", 6))
crystal_shape = input_data.get("crystal_shape", "fcc")

# Get atoms per unit cell
if GB:
    try:
        total_atoms = len(read("POSCAR-gb"))
    except Exception as e:
        raise FileNotFoundError(
            f"!! Grain boundary is enabled (`GB=True`), but POSCAR-gb was not found.\n"
            "Please provide a valid POSCAR-gb file in the working directory. !!"
        )
else:
    try:
        atoms = bulk('X', crystalstructure=crystal_shape, cubic=True)  # Ensuring orthorhombic setting
        atoms_per_unit_cell = len(atoms)
    except Exception as e:
        atoms = bulk('X', crystalstructure=crystal_shape, orthorhombic=True)  # Ensuring orthorhombic setting
        atoms_per_unit_cell = len(atoms)

    

    total_atoms = int(size**3 * atoms_per_unit_cell)

# Parse additives
undoped = False
if input_data["additives"].lower() == "none":
    additives = []
    undoped = True
else:
    cleaned_value = input_data["additives"].strip()[1:-1]
    additives = [(x.strip() for x in add.strip('()').split(',')) for add in cleaned_value.split('), (')]

# Function to update the input file
def update_input_file(file_path, updates):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for key, value in updates.items():
        for i, line in enumerate(lines):
            if line.startswith(key):
                if isinstance(value, dict):
                    value_str = ', '.join([f"{k}={v}" for k, v in value.items()])
                elif isinstance(value, list):
                    value_str = "[" + ", ".join(f"({x[0]}, {x[1]}, {x[2]})" for x in value) + "]"
                else:
                    value_str = str(value)
                lines[i] = f"{key}: {value_str}\n"
                break

    with open(file_path, 'w') as file:
        file.writelines(lines)

# Function to execute a single simulation
def run_simulation(fraction, run_id=None):
    num_dopants = round(fraction * total_atoms)
    updated_additives = [(host, dopant, num_dopants) for host, dopant, _ in additives]

    # Generate a unique directory name
    if undoped:
        base_dir_name = "Alloy"
    else:
        base_dir_name = f"Alloy_{'_'.join([dopant for _, dopant, _ in updated_additives])}_{int(fraction * 100)}"

    # Add a counter if directory already exists
    dir_name = base_dir_name
    counter = 1
    while os.path.exists(dir_name):
        dir_name = f"{base_dir_name}_run{counter}"
        counter += 1

    os.makedirs(dir_name, exist_ok=True)

    # Copy relevant files to the new directory
    items_to_copy = ["chgnet_src", "src", "hybrid_mcmd.py", "input_file"]
    for item in items_to_copy:
        src_path = os.path.abspath(item)
        dest_path = os.path.join(dir_name, os.path.basename(item))
        if os.path.isdir(src_path):
            shutil.copytree(src_path, dest_path, dirs_exist_ok=True)
        elif os.path.isfile(src_path):
            shutil.copy(src_path, dest_path)

    # Change to the unique directory
    os.chdir(dir_name)

    # Update input file for the specific dopant fraction
    updates = {
        "composition": {k.strip(): float(v.strip()) for k, v in [item.split('=') for item in input_data['composition'].split(',')]},
        "grain_boundary": GB,
        "additives": updated_additives if not undoped else "None"
    }

    update_input_file("input_file", updates)

    # Run the simulation
    os.system("python3 hybrid_mcmd.py")

    # Return to the parent directory
    os.chdir("..")

    print(f"Simulation for fraction {fraction * 100}% completed in {dir_name}.")

# Main execution with parallelization
if __name__ == "__main__":
    if batch_fractions:
        with Pool(processes=cpu_count()) as pool:
            pool.map(run_simulation, batch_fractions)

    else:
        print("Please use serial version of execute.py for single simulations.")
