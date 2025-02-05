#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:25:49 2024

@author: dolez
"""

import os
import shutil

# FCC supercell with 6x6x6x4 = 864 atoms
# BCC supercell with 7x7x7x2 = 686 atoms

atom_count = {"CrNi" : 864, 
              "CrFe" : 686}

# Path to your execute.py script
script_path = "hybrid_mc_md.py"

# Files and directories to copy
files_to_copy = ["POSCAR-1", "data", "structures"]

system = "CrNi"

# Metal host composition
total_atoms = atom_count[system]  

undoped = False
# if grain boundary structure (True)
GB = True

# Loop through atomic fractions of dopants
dopant_fractions = [0.20, 0.10, 0.08, 0.06, 0.04, 0.02, 0.01]
metal_host = 'Cr'  # Example metal host
dopant = 'B'

# label the runs 1-5
if undoped or GB:
    dopant_fractions = [0.01, 0.02, 0.03, 0.04, 0.05]

for fraction in dopant_fractions:
    num_dopants = round(fraction * total_atoms)  # Calculate number of dopants
    
    if GB and not undoped:
        num_dopants = 50
    
    # Read the original script
    with open(script_path, 'r') as file:
        lines = file.readlines()
    
    # Update the additives line
    for i, line in enumerate(lines):
        if line.startswith("additives ="):
            lines[i] = f"additives = ('{metal_host}', '{dopant}', {num_dopants})\n"
            
            if undoped:
                lines[i] = "additives = None\n"
            
        
        if line.startswith("grain_boundary ="):
            lines[i] = f"grain_boundary = {GB}\n"
    
    with open(script_path, 'w') as file:
        file.writelines(lines)
    
    file.close()
    
    # Execute the script
    os.system(f"python3 {script_path}")
    
    # Record the simulation data
    
    if undoped:
        # Create directory for this configuration
        dir_name = f"{system}_{int(fraction * 100)}"
        os.makedirs(dir_name, exist_ok=True)
    
    else:
        # Create directory for this configuration
        dir_name = f"{system}+{dopant}_{int(fraction * 100)}"
        os.makedirs(dir_name, exist_ok=True)
    
    
    
    # Copy required files and directories into the directory
    for item in files_to_copy:
        if os.path.exists(item):  # Check if the file or directory exists
            if os.path.isdir(item):  # Copy directory
                shutil.copytree(item, os.path.join(dir_name, os.path.basename(item)), dirs_exist_ok=True)
                shutil.rmtree(item) # remove directory so new one will be created
            else:  # Copy file
                shutil.copy2(item, os.path.join(dir_name, os.path.basename(item)))
        else:
            print(f"Warning: {item} not found. Skipping.")
        
        