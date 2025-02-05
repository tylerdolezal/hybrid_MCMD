from ase.io import read, write
import os
import numpy as np
import pandas as pd
import random
import copy
import json
import src.myfuncs as fun
import src.lammps_functions as lfun
import src.ht_sampling_routine as ht

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
            - O2: bool, whether to include oxygen at half the dopant concentration
            - vacancies: bool, whether to include vacancies

    Returns:
        None
    """
    # Validate input
    required_keys = ['composition', 'crystal_shape', 'grain_boundary', 'md_params',
                     'num_mc_steps', 'md_interval', 'size', 'supcomp_command',
                     'continue_run', 'additives', 'O2', 'vacancies']

    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required configuration key: {key}")

    # Create required directories
    for directory in ["structures", "data", "data/epochs"]:
        os.makedirs(directory, exist_ok=True)

    _, temperature, potential_type = config['md_params']
    lfun.update_md_input_file(config['md_params'])

    if config['md_params'][2] == "chgnet":
        lfun.update_chgnet_input_file(config['md_params'])
    else:
        lfun.update_md_input_file(config['md_params'])

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
            config['O2'],
            config['crystal_shape'],
            config['size'],
            config['vacancies']
        )

        # Initial energy setup
        E0 = np.loadtxt("energy.txt")
        energies = [E0, E0]
        np.savetxt('data/energies', energies)

        # Initialize move data DataFrame
        columns = ['Atom Type', 'r_0', 'r_f', 'Times Moved']
        PFM_data = pd.DataFrame(columns=columns)

    else:
        system = read('POSCAR-1', format='vasp')
        PFM_data = pd.read_csv('data/phase_field_model_data.csv', index_col=None)
        lfun.update_modfiles(config['md_params'])

    if config['continue_run']:
        energies = list(np.loadtxt('data/energies'))

        # Parse MC statistics if in continue run mode
        save, accepted, rejected, move_counts, steps_completed = fun.parse_mc_statistics()

        with open('data/species_counter.json', 'r') as f:
            species_counter = json.load(f)

    # Main simulation logic
    move_list = ['swap']
    if config['additives']:
        move_list = ['swap', 'new_host', 'new_host', 'shuffle']

    snapshot_every = 100
    num_mc_steps = config['num_mc_steps']
    md_interval = config['md_interval']

    for mc_step in range(steps_completed + 1, num_mc_steps + 1):

        # Randomly choose a move type
        move_type = random.choice(move_list)

        # Select atoms for the move
        if move_type not in ['cluster', 'shuffle', 'cluster_hop']:
            swap_pairs = fun.select_random_atoms(system, move_type)
        else:
            swap_pairs = (0, 0)  # Placeholder for unsupported moves

        if save >= md_interval:
            print(f'Running HT MC sampling routine at step {mc_step}!')
            move_type = 'MD'
            delta_E, new_energy, swap_pairs = ht.high_temperature_sampling(
                system, energies[-1], config['md_params'], config['additives']
            )

            # Automatically accept the new state
            old_system = copy.deepcopy(system)
            system = read('CONTCAR', format='vasp')

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

            with open('data/species_counter.json', 'w') as file:
                json.dump(species_counter, file)

        else:
            delta_E, new_energy = fun.calculate_energy_change(
                system, energies, swap_pairs, move_type, False, config['supcomp_command']
            )

            if random.random() < np.exp(-delta_E / (k_B * temperature)):
                old_system = copy.deepcopy(system)
                system = read('CONTCAR', format='vasp')

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

                with open('data/species_counter.json', 'w') as file:
                    json.dump(species_counter, file)

            else:
                energies.append(energies[-1])
                rejected += 1

        with open('data/MonteCarloStatistics', 'w') as file:
            file.write(f'Steps Completed: {mc_step}\n')
            file.write(f'Acceptance %: {accepted/mc_step * 100}\n')
            file.write(f'Rejection %: {rejected/mc_step * 100}\n')
            file.write(f"\nAccepted Swaps: {move_counts['swap']}\n")
            file.write(f"New Hosts Accepted: {move_counts['new_host']}\n")
            file.write(f"Translates Accepted: {move_counts['translate']}\n")
            file.write(f"Cluster Hops Accepted: {move_counts['cluster_hop']}\n")
            file.write(f"Cluster Shuffles Accepted: {move_counts['shuffle']}\n")
            file.write(f"MD Simulations Accepted: {move_counts['MD']}\n")

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
            key, value = line.strip().split('=')
            config[key.strip()] = eval(value.strip())
    return config

# To run the routine
example_config = read_config_file()
hybrid_md_mc_routine(example_config)