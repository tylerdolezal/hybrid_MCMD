import os
import shutil
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
if input_data["additives"].lower() == "none":
    additives = []
    undoped = True
else:
    cleaned_value = input_data["additives"].strip()[1:-1]
    additives = [(x.strip() for x in add.strip('()').split(',')) for add in cleaned_value.split('), (')]

# If executing a batch job
if batch_fractions:
    print(f"Running batch execution for fractions: {batch_fractions}")

    for fraction in batch_fractions:
        num_dopants = round(fraction * total_atoms)
        updated_additives = [(host, dopant, num_dopants) for host, dopant, _ in additives]

        # Update input file for the current fraction
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

        updates = {
            "composition": {k.strip(): float(v.strip()) for k, v in [item.split('=') for item in input_data['composition'].split(',')]},
            "grain_boundary": GB,
            "additives": updated_additives if not undoped else "None"
        }

        update_input_file("input_file", updates)

        # Run the simulation
        os.system("python3 hybrid_mcmd.py")

        # Organize results
        # Base directory name
        base_dir_name = "Alloy"+f"{'_'.join([dopant for _, dopant, _ in updated_additives])}{int(fraction * 100)}"

        # Check if the directory already exists, and add a counter if necessary
        dir_name = base_dir_name
        counter = 1
        while os.path.exists(dir_name):
            dir_name = f"{base_dir_name}_run{counter}"
            counter += 1

        # Create the final unique directory
        os.makedirs(dir_name, exist_ok=True)

        # Copy result files
        for item in ["POSCAR-1", "data", "structures"]:
            if os.path.exists(item):
                dest_path = os.path.join(dir_name, item)
                if os.path.isdir(item):
                    shutil.copytree(item, dest_path)
                else:
                    shutil.copy(item, dest_path)

        print(f"Simulation for {', '.join([dopant for _, dopant, _ in updated_additives])} fraction {fraction*100}% completed.")

else:
    # If no batch job, execute normally
    print("Executing standard single-run hybrid MCMD...")
    os.system("python3 hybrid_mcmd.py")
