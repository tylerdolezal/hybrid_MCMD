#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 15:56:19 2024

@author: dolez
"""
import copy
import src.myfuncs as fun
import src.lammps_functions as lfun

def hybrid_md_sampling(system, current_energy, md_params, additives):
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
    # make sure we swap for chgnet runs
    if md_params[2] == 'chgnet':
        import chgnet_src.myfuncs as fun

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

