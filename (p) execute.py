import os
import shutil
from multiprocessing import Pool, cpu_count

# Configuration settings
system = "CrNi"
undoped = False

# Read and parse input file
with open("input_file", 'r') as file:
    input_data = {line.split(':')[0].strip(): line.split(':')[1].strip() for line in file if ':' in line}

GB = input_data["grain_boundary"].lower() == 'true'
additives = input_data['additives'].split(',')
metal_host = additives[0].strip()
dopant = additives[1].strip()

if input_data["additives"] == "None":
    undoped = True

size = int(input_data['size'])
crystal_shape = input_data['crystal_shape']
cuc = {'fcc': 4, 'bcc': 2}
total_atoms = int(size**3 * cuc[crystal_shape])

if GB and not os.path.exists("POSCAR-gb"):
    raise FileNotFoundError("POSCAR-gb file is required but not found. Please ensure the file exists in the working directory.")

dopant_fractions = [0.20, 0.10, 0.08, 0.06, 0.04, 0.02, 0.01]
if undoped or GB:
    dopant_fractions = [0.01, 0.02, 0.03, 0.04, 0.05]

# Function to update the input file
def update_input_file(file_path, updates):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for key, value in updates.items():
        for i, line in enumerate(lines):
            if line.startswith(key):
                if isinstance(value, dict):
                    value_str = ', '.join([f"{k}={v}" for k, v in value.items()])
                    lines[i] = f"{key}: {value_str}\n"
                elif isinstance(value, tuple):
                    value_str = ', '.join(map(str, value))
                    lines[i] = f"{key}: {value_str}\n"
                else:
                    lines[i] = f"{key}: {value}\n"
                break

    with open(file_path, 'w') as file:
        file.writelines(lines)

# Function to execute a single simulation
def run_simulation(fraction):
    # Unique directory name for each run
    run_dir = f"run_{int(fraction * 100)}"
    os.makedirs(run_dir, exist_ok=True)

    # Files and directories to copy into run_dir
    items_to_copy = ["chgnet_src", "src", "hybrid_mcmd.py", "input_file"]

    # Copy specified files and directories into the run_dir
    for item in items_to_copy:
        src_path = os.path.abspath(item)
        dest_path = os.path.join(run_dir, os.path.basename(item))
        if os.path.isdir(src_path):
            shutil.copytree(src_path, dest_path, dirs_exist_ok=True)
        elif os.path.isfile(src_path):
            shutil.copy(src_path, dest_path)

    # Change to the unique directory
    os.chdir(run_dir)

    # Update input_file with specific dopant fraction
    num_dopants = round(fraction * total_atoms)
    updates = {
        "composition": {k.strip(): float(v.strip()) for k, v in [item.split('=') for item in input_data['composition'].split(',')]},
        "grain_boundary": GB,
        "additives": f"{metal_host}, {dopant}, {num_dopants}"
    }
    if undoped:
        updates["additives"] = None

    update_input_file("input_file", updates)

    # Run the simulation
    os.system("python3 hybrid_mcmd.py")

    # Return to the parent directory
    os.chdir("..")

    print(f"Simulation for {dopant} fraction {fraction * 100}% completed in {run_dir}.")

# Main execution with parallelization
if __name__ == "__main__":
    with Pool(processes=cpu_count()) as pool:
        pool.map(run_simulation, dopant_fractions)
