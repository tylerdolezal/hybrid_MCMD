from ase.io import read, write
import os
import numpy as np
import pandas as pd
import random
import copy
import json
import src.lammps_functions as lfun
import src.hybrid_md_routine as ht

# Constants
k_B = 8.617333262145e-5  # Boltzmann constant in eV/K

def hybrid_md_mc_routine(config):
    """
    Executes the hybrid MD-MC routine using parameters provided in a configuration dictionary.

    Parameters:
        config (dict): Configuration dictionary containing the following keys:
            - composition: dict, element composition of the system (e.g., {'Ni': 0.7, 'Cr': 0.3})
            - crystal_shape: str, shape of the crystal (e.g., 'fcc', 'bcc')
            - grain_boundary: bool or dict, grain boundary parameters (True for default, or provide a dictionary with specific parameters)
            - md_params: tuple, MD simulation parameters (e.g., (0, 1073, 'eam'))
            - num_mc_steps: int, number of Monte Carlo steps
            - md_interval: int, MC steps between MD steps
            - size: int or tuple, dimensions of the simulation box (e.g., 6 for cubic or (10, 10, 10))
            - supcomp_command: str, supercomputer command to run simulations
            - continue_run: bool, whether to continue from a previous run
            - additives: tuple, dopant type and count (e.g., ('B', 10))
            - vacancies: bool, whether to include vacancies

    Returns:
        None
    """
    import src.myfuncs as fun
    # Validate input
    required_keys = ['composition', 'crystal_shape', 'grain_boundary', 'randomize', 'md_params',
                     'num_mc_steps', 'md_interval', 'size', 'supcomp_command',
                     'continue_run', 'additives', 'vacancies', 'metal_library', 'surface', 'local_swap']

    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required configuration key: {key}")

    local = config['local_swap']
    # Create required directories
    for directory in ["structures", "data", "data/epochs"]:
        os.makedirs(directory, exist_ok=True)

    _, temperature, potential_type = config['md_params']
    lfun.update_md_input_file(config['md_params'])

    # set global freeze threshold based on the surface input
    if config['surface']:
        threshold = fun.set_global_threshold(config['surface'])
        lfun.update_md_input_file(config['md_params'], threshold)
        
    # update our list of interstitials based on the additives
    # and surface adsorbates
    if config['additives']:
        species = []
        for dopant in config['additives']:
            species.append(dopant[1])
        fun.set_interstitials(species)

    if config['surface']:
        species = [config['surface'][0]]
        fun.set_interstitials(species, surface=True)        
    
    # switch to chgnet myfuncs if using chgnet
    if potential_type == 'chgnet':
        import chgnet_src.myfuncs as fun

    if not config['continue_run']:
        # Check for POSCAR-gb if GB is True
        if config['grain_boundary'] and not os.path.exists("POSCAR-gb"):
            raise FileNotFoundError("POSCAR-gb file not found for grain boundary configuration.")

        system = fun.initialize_system(
            config['composition'],
            config['grain_boundary'],
            config['supcomp_command'],
            config['md_params'],
            config['additives'],
            config['crystal_shape'],
            config['size'],
            config['vacancies'],
            config['randomize'],
            config['surface']
        )

        # Initial energy setup
        E0 = np.loadtxt("energy.txt")
        energies = [E0, E0]
        np.savetxt('data/energies', energies)

        # Initialize move data DataFrame
        columns = ['Atom Type', 'r_0', 'r_f', 'Times Moved']
        PFM_data = pd.DataFrame(columns=columns)

        # initialize our files to avoid errors in restarting
        species_counter = fun.initialize_mc_step_counter(system)
        with open('data/species_counter.json', 'w') as file:
            json.dump(species_counter, file)
        PFM_data.to_csv('data/phase_field_model_data.csv', index=False)

        accepted = 0; rejected = 0
        steps_completed = 0
        save = 0
        move_counts = {'swap': 0, 'new_host': 0, 'flip': 0, 'shuffle': 0, 'MD': 0}

    else:
        system = read('POSCAR-1', format='vasp')
        PFM_data = pd.read_csv('data/phase_field_model_data.csv', index_col=None)

        energies = list(np.loadtxt('data/energies'))

        # Parse MC statistics if in continue run mode
        save, accepted, rejected, move_counts, steps_completed = fun.parse_mc_statistics()

        with open('data/species_counter.json', 'r') as f:
            species_counter = json.load(f)


    write("POSCAR", system, format='vasp', direct=True, sort=False)
    lfun.update_modfiles(config['md_params'])

    # Main simulation logic
    move_list = ['swap']
    metal_choices = config['metal_library']
    if config['additives'] or config['surface']:
        move_list = ['swap', 'new_host', 'new_host', 'shuffle']
        if metal_choices:
            move_list =  move_list + ['flip']

    snapshot_every = 100
    num_mc_steps = config['num_mc_steps']
    md_interval = config['md_interval']
    update_list = False
    for mc_step in range(steps_completed + 1, num_mc_steps + 1):

        # One-time move_list update for surface calculations; give time for dimer 
        # dissociation before attempting interstitial swaps at surface
        if not update_list and mc_step > 1000 and config['surface']:
            move_list = fun.update_move_list(config['surface'])

            if metal_choices:
                move_list += ['flip']

            update_list = True

        # Randomly choose a move type
        move_type = random.choice(move_list)

        # Select atoms for the move
        if move_type not in ['shuffle', 'flip']:
            swap_pairs = fun.select_random_atoms(system, move_type, local)
        else:
            swap_pairs = (0, 0)  # Placeholder for unsupported moves

        if save >= md_interval:
            print(f'Running HT MC sampling routine at step {mc_step}!')
            move_type = 'MD'
            delta_E, new_energy, swap_pairs = ht.hybrid_md_sampling(
                system, energies[-1], config['md_params'], config['additives']
            )

            # Automatically accept the new state
            old_system = copy.deepcopy(system)
            system = read('CONTCAR', format='vasp')

            PFM_data, species_counter = fun.update_phase_field_dataset(
                PFM_data, old_system, system, delta_E, swap_pairs, species_counter, move_type
            )

            write("POSCAR-1", system, format='vasp', direct=True, sort=False)
            move_counts[move_type] += 1

            if save >= md_interval:
                save = 0

            energies.append(new_energy)
            accepted += 1
            save += 1

            with open('data/species_counter.json', 'w') as file:
                json.dump(species_counter, file)

        else:
            delta_E, new_energy = fun.calculate_energy_change(
                system, energies, swap_pairs, move_type, False, config['supcomp_command'], local, metal_choices
            )

            if random.random() < np.exp(-delta_E / (k_B * temperature)):
                old_system = copy.deepcopy(system)
                system = read('CONTCAR', format='vasp')

                if move_type != 'flip':
                    PFM_data, species_counter = fun.update_phase_field_dataset(
                        PFM_data, old_system, system, delta_E, swap_pairs, species_counter, move_type
                    )

                write("POSCAR-1", system, format='vasp', direct=True, sort=True)
                move_counts[move_type] += 1

                if save >= md_interval:
                    save = 0

                energies.append(new_energy)
                accepted += 1
                save += 1

                species_counter = fun.add_new_species(system, species_counter)

                with open('data/species_counter.json', 'w') as file:
                    json.dump(species_counter, file)

            else:
                energies.append(energies[-1])
                rejected += 1

        with open('data/MonteCarloStatistics', 'w') as file:
            file.write(f'Steps Completed: {mc_step}\n')
            file.write(f'Acceptance %: {accepted/mc_step * 100}\n')
            file.write(f'Rejection %: {rejected/mc_step * 100}\n')
            file.write(f"\nAccepted Metal Swaps: {move_counts['swap']}\n")
            file.write(f"Accepted Int Swaps: {move_counts['swap_ints']}\n")
            file.write(f"New Hosts Accepted: {move_counts['new_host']}\n")
            file.write(f"Flips Accepted: {move_counts['flip']}\n")
            file.write(f"Cluster Shuffles Accepted: {move_counts['shuffle']}\n")
            file.write(f"MD Simulations Accepted: {move_counts['MD']}\n")
            file.write(f'Steps for MD: {save}')

        file.close()

        fun.snapshots(mc_step, snapshot_every)

        np.savetxt("data/energies", energies)

    print("Simulation complete.")

# Function to read configuration from mcmd.in file
def read_config_file():
    """
    Reads the hardcoded mcmd.in configuration file in a text format and parses it into a dictionary.

    Returns:
        dict: Configuration dictionary
    """
    config = {}
    with open('input_file', 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()

                if key == 'composition':
                    config[key] = {item.split('=')[0]: float(item.split('=')[1]) for item in value.split(', ')}
                elif key == 'md_params':
                    config[key] = tuple(map(lambda x: int(x) if x.isdigit() else x, value.split(', ')))
                elif key in {'num_mc_steps', 'md_interval', 'size'}:
                    config[key] = int(value)
                elif key in {'grain_boundary', 'continue_run', 'vacancies', 'randomize', 'local_swap'}:
                    config[key] = value.lower() == 'true'
                elif key == 'additives':
                    if value.strip().lower() == "none":
                        config[key] = None
                    
                    else:
                        # Remove outer brackets and strip spaces
                        cleaned_value = value.strip().strip("[]")

                        # Split into individual tuples based on "), ("
                        tuple_strings = cleaned_value.split("), (")

                        additives_list = []
                        for tuple_str in tuple_strings:
                            # Remove any stray parentheses and split by commas
                            elements = tuple_str.replace("(", "").replace(")", "").split(",")
                            
                            # Ensure we have exactly 3 elements
                            if len(elements) != 3:
                                raise ValueError(f"Unexpected format for additives: {value}")
                            
                            # Convert the third element to an integer
                            additive_tuple = (elements[0].strip(), elements[1].strip(), int(elements[2].strip()))
                            
                            # Append to list
                            additives_list.append(additive_tuple)

                        # Store as list of tuples
                        config[key] = additives_list

                elif key == 'metal_library':
                    value = value.strip()

                    # If explicitly "None", store as None
                    if value.lower() == "none":
                        config[key] = None
                    else:
                        # Remove brackets and strip spaces
                        cleaned_value = value.replace("[", "").replace("]", "").strip()
                        
                        # Split by commas, strip spaces, remove extra quotes, and filter out empty entries
                        config[key] = [v.strip().strip("'").strip('"') for v in cleaned_value.split(',') if v.strip()]

                elif key == 'surface':
                    if value.strip().lower() == "none":
                        config[key] = None
                    else:
                        # Remove square brackets and strip whitespace
                        cleaned_value = value.replace("[", "").replace("]", "").strip()
                        
                        # Split by commas, clean up entries
                        parsed_values = [v.strip().strip("'").strip('"') for v in cleaned_value.split(',') if v.strip()]

                        # Convert last element to float or int if possible
                        if parsed_values and parsed_values[-1].replace('.', '', 1).isdigit():
                            parsed_values[-1] = float(parsed_values[-1])
                        
                        config[key] = parsed_values
                else:
                    config[key] = value
    return config


# To run the routine
example_config = read_config_file()
hybrid_md_mc_routine(example_config)