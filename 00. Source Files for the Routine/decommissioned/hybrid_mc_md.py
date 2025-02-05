#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:09:24 2024

@author: dolez
"""
import os
import json
import copy
import random
import numpy as np
import pandas as pd
import src.myfuncs as fun
from ase.io import read, write
import src.ht_sampling_routine as ht
import src.lammps_functions as lfun

# Constants
k_B = 8.617333262145e-5  # Boltzmann constant in eV/K

# Hybrid MD-MC Routine
def hybrid_md_mc_routine(composition, crystal_shape, grain_boundary, md_params, num_mc_steps, md_interval, size, supcomp_command, continue_run, additives, O2, vacancies):

    if not os.path.exists("structures"):
        os.makedirs("structures")

    if not os.path.exists("data"):
        os.makedirs("data")

    if not os.path.exists("data/epochs"):
        os.makedirs("data/epochs")

    temperature = md_params[1]

    lfun.update_md_input_file(md_params)

    if not continue_run:
        system = fun.initialize_system(composition, grain_boundary, supcomp_command, md_params, additives, O2, crystal_shape, size, vacancies)

        # initial cell is relaxed within initialize_system function
        E0 = np.loadtxt("energy.txt")
        # ensure the array has a len for grow cluster check
        energies = [E0,E0]

        np.savetxt('data/energies', energies)

        # Initialize a DataFrame to store the move data
        columns = ['Atom Type', 'r_0', 'r_f', 'Times Moved']
        PFM_data = pd.DataFrame(columns=columns)

    else:
        system = read('POSCAR-1', format='vasp')

        PFM_data = pd.read_csv('data/phase_field_model_data.csv', index_col=None)

        lfun.update_modfiles(md_params)





    accepted = 0; rejected = 0
    steps_completed = 0
    save = 0
    move_counts = {'swap': 0, 'new_host': 0, 'translate': 0, 'cluster_hop': 0, 'shuffle': 0, 'MD': 0}

    species_counter = fun.initialize_mc_step_counter(system)

    if continue_run:
        energies = list(np.loadtxt('data/energies'))

        with open('data/species_counter.json', 'r') as file:
            species_counter = json.load(file)


        save, accepted, rejected, move_counts, steps_completed = fun.parse_mc_statistics()

    move_list = ['swap']
    if additives:
        #move_list = ['swap', 'new_host', 'translate', 'cluster', 'shuffle']
        # equal chances to move a metal or interstitial
        move_list = ['swap', 'new_host', 'new_host', 'shuffle']

    snapshot_every = 100
    for mc_step in range(steps_completed+1, num_mc_steps+1):

        # Randomly choose between swap and translational move
        move_type = random.choice(move_list)
        # Grab the indices of metals who get swapped/translated
        # or the metal, B/C atom for new host pair
        if move_type not in ['cluster', 'shuffle', 'cluster_hop']:
            swap_pairs = fun.select_random_atoms(system, move_type)

        else:
            swap_pairs = (0,0) # not used for cluster translations and shuffle


        if save >= md_interval:
            np.savetxt(f'Running HT MC sampling routine {mc_step}!', [])
            # set move type to MD so all we are trying is new thermal configuration
            move_type = 'MD'
            delta_E, new_energy, swap_pairs = ht.high_temperature_sampling(system, energies[-1], md_params, additives)

        else:
            delta_E, new_energy = fun.calculate_energy_change(system, energies, swap_pairs, move_type, False, supcomp_command)

        if random.random() < np.exp(-delta_E / (k_B * temperature)):

            old_system = copy.deepcopy(system)

            # the CONTCAR file was last created for calculating new_energy,
            # so by reading in this file we are updated system to the accepted state
            system = read('CONTCAR', format='vasp')

            # grab the final configutaion details
            PFM_data, species_counter = fun.update_phase_field_dataset(PFM_data, old_system, system, delta_E, swap_pairs, species_counter, move_type)

            # for book keeping dump the latest accepted state in POSCAR-1
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

            if save >= md_interval:
                save = 0

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
            file.write(f'Steps for MD: {save}')
        file.close()

        fun.snapshots(mc_step, snapshot_every)

        np.savetxt("data/energies", energies)

        if mc_step % snapshot_every == 0:
                os.system(f'cp POSCAR-1 structures/POSCAR-{mc_step}')




# Parameters
num_mc_steps = 6000  # Total number of MC steps
md_interval = 6001  # Number of successful swaps before running MD
md_params = (0, 1073, 'eam')  # number of md steps and Temp in K
# these are turned down for now since I'm starting with no MD

composition = {'Ni': 0.58, 'Cr': 0.22, 'Fe': 0.05, 'Mo': 0.10, 'Nb': 0.03, 'Ti': 0.02}


composition = {'Ni': 0.70, 'Cr': 0.30}

crystal_shape = 'fcc'
supercell_multiplier = 6

# metal host, dopant, number of dopants to add
additives = ('B', 10)
# whether or not to through 1/2 the amount of O into the system
O2 = False
vacancies = False
continue_run = False
grain_boundary = True

# command to execute LAMMPS on the system
supcomp_command = "mpiexec -n 128 mylammps/build/lmp" # for warhawk

# supcomp_command = "lmp_serial" # for Matlantis

# Run the hybrid MD-MC routine
hybrid_md_mc_routine(composition, crystal_shape, grain_boundary, md_params, num_mc_steps, md_interval, supercell_multiplier, supcomp_command, continue_run, additives, O2, vacancies)
