#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 15:56:19 2024

@author: dolez
"""
import os
import copy
import random
import numpy as np
from ase.io import read, write
import src.myfuncs as fun
import src.lammps_functions as lfun

# Constants
k_B = 8.617333262145e-5  # Boltzmann constant in eV/K

def high_temperature_sampling(system, current_energy, md_params, additives):
    """
    Perform high temperature sampling by running an MD simulation, saving the post-MD structure,
    relaxing that structure, and returning the energy difference.

    Parameters:
        system (ASE Atoms object): The atomic configuration for simulation.
        current_energy (float): The current energy of the system.
        md_params (tuple): MD simulation parameters.
        additives (tuple): Dopant type and count.

    Returns:
        delta_E (float): The energy difference after relaxation.
        new_energy (float): The new energy after relaxation.
        swap_pairs (list): The swap pairs used in the simulation.
    """
    # Start by heating the current system up
    system_copy = copy.deepcopy(system)
    system_copy, E0 = fun.run_md_simulation(system_copy, md_params[2])

    # Relax the new structure
    lfun.update_md_input_file(md_params)
    system_copy, Ef = fun.relax(system_copy, 'new_host', md_params[2])

    # Calculate the energy difference
    delta_E = Ef - current_energy

    # Return the energy difference, new energy, and an empty swap_pairs list
    return delta_E, Ef, []

# Commented out the old high temperature sampling routine
"""
def high_temperature_sampling(system, current_energy, md_params, additives):
    
    ht_move_counts = {'swap': 0, 'new_host': 0, 'translate': 0, 'cluster_hop': 0, 'shuffle': 0, 'MD': 0}
    
    accepted, rejected = 0, 0
    
    temperature = md_params[1]
    os.system("cp relax-ht.in relax.in")
    system_copy = copy.deepcopy(system)
    
    # start by heating the current system up
    system_copy, E0 = fun.run_md_simulation(system_copy)
    # now we have the HT structure, with HT fluctuations present
    
    # attempt 100 new moves on the HT structure
    move_list = ['swap']
    if additives:
        move_list = ['swap', 'new_host', 'shuffle']
        
    ht_energies = [E0]
    success = 0
    for ht_step in range(1,101):
            
        # Randomly choose between swap and translational move
        move_type = random.choice(move_list)
        # Grab the indices of metals who get swapped/translated
        # or the metal, B/C atom for new host pair
        if move_type not in ['cluster', 'shuffle', 'cluster_hop']:
            swap_pairs = fun.select_random_atoms(system_copy, move_type)
        
        else:
            swap_pairs = (0,0) # not used for cluster translations and shuffle

        delta_E, new_energy = fun.calculate_energy_change(system_copy, ht_energies, swap_pairs, move_type, False)
            
        if random.random() < np.exp(-delta_E / (k_B * temperature)):
                        
            system_copy = read('CONTCAR', format='vasp')
            
            ht_energies.append(new_energy)
            
            success += 1
            accepted += 1
            ht_move_counts[move_type] += 1
        
        else:
            ht_energies.append(ht_energies[-1])
            rejected += 1
            
        np.savetxt('success_at_ht', [success])
        
        with open('data/Monte Carlo Statistics (HT)', 'w') as file:
            file.write(f'Steps Completed: {ht_step}\n')
            file.write(f'Acceptance %: {accepted/ht_step * 100}\n')
            file.write(f'Rejection %: {rejected/ht_step * 100}\n')
            file.write(f"\nAccepted Swaps: {ht_move_counts['swap']}\n")
            file.write(f"New Hosts Accepted: {ht_move_counts['new_host']}\n")
            file.write(f"Translates Accepted: {ht_move_counts['translate']}\n")
            file.write(f"Cluster Hops Accepted: {ht_move_counts['cluster_hop']}\n")
            file.write(f"Cluster Shuffles Accepted: {ht_move_counts['shuffle']}\n")
            file.write(f"MD Simulations Accepted: {ht_move_counts['MD']}\n")
        file.close()
    
    # cool the final structure down
    # revert back to CG style relax.in file
    lfun.update_md_input_file(md_params)
    system_copy, Ef = fun.relax(system_copy)
    
    # make sure we give back proper swap_pairs for the
    # phase model data set
    if move_type == 'shuffle':
        with open('shuffle_pairs') as file:
            line = file.readline().strip()  # Read the first line and strip any extra whitespace/newline characters
            parts = line.split(',')  # Split the line at the comma to separate the indices
            swap_pairs = [tuple(map(int, parts))]
    
    # give back the proper cluster hop pairs for HT MD attempts
    if move_type == 'cluster_hop':
        with open('cluster_hop_pairs') as file:
            line = file.readline().strip()  # Read the first line and strip any extra whitespace/newline characters
            parts = line.split(',')  # Split the line at the comma to separate the indices
            swap_pairs = [tuple(map(int, parts))]
    
    return Ef - current_energy, Ef, swap_pairs
"""