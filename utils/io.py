import os
import numpy as np
import pandas as pd
from ase.io import read, write

def save_snapshot(sim):
    """
    Handles periodic checkpoints, saving structures and updating statistics.
    """
    # 1. Save current POSCAR for restart capability 
    write("POSCAR-1", sim.atoms, format='vasp', direct=True, sort=True)
    
    # 2. Update energy log 
    np.savetxt("data/energies", sim.energies)

    # 3. Create a backup of the current structure 
    os.system(f'cp POSCAR-1 structures/POSCAR-{sim.mc_step}')
    
    # 4. Write MonteCarloStatistics in the required text format 
    write_mc_statistics(sim)

def record_spectral_entry(energy, dEf, dEr, system, n_solute, solute, output_path="data/spectral_log.csv"):
    """
    Records detailed metadata for spectral/NEB mapping.
    
    Parameters:
    - energy: Total potential energy of the configuration.
    - dEf/dEr: Forward and reverse NEB activation barriers.
    - system: The ASE Atoms object containing the decorated structure.
    - n_solute: The current number of solute neighbors.
    - solute: The symbol of the solute element.
    """
    # Identify the interstitial atom (the last atom added/moved)
    # We assume the last index is our target dopant
    idx = len(system) - 1
    pos = system[idx].position
    
    # Build the row dictionary
    row = {
        "energy (eV)": energy,
        "dE_forward (eV)": dEf,
        "dE_reverse (eV)": dEr,
        "x (Å)": pos[0],
        "y (Å)": pos[1],
        "z (Å)": pos[2],
        f"n_{solute}": n_solute
    }

    # Use a DataFrame to handle CSV appending
    df = pd.DataFrame([row])
    header = not os.path.exists(output_path)
    df.to_csv(output_path, mode='a', index=False, header=header)

def write_mc_statistics(sim):
    """
    Writes statistics to data/MonteCarloStatistics. 
    Maintains exact string formatting for compatibility with parse_mc_statistics.
    """
    path = 'data/MonteCarloStatistics'
    with open(path, 'w') as f:
        f.write(f'Steps Completed: {sim.mc_step}\n')
        f.write(f'Acceptance %: {(sim.accepted / sim.mc_step) * 100 if sim.mc_step > 0 else 0.0}\n')
        f.write(f'Rejection %: {(sim.rejected / sim.mc_step) * 100 if sim.mc_step > 0 else 0.0}\n\n')
        
        f.write(f'Accepted Moves\n')
        f.write(f"Metal Swaps: {sim.move_stats['accepted']['swap']}\n")
        f.write(f"Interstitial Swaps: {sim.move_stats['accepted']['swap_ints']}\n")
        f.write(f"New Hosts: {sim.move_stats['accepted']['new_host']}\n")
        f.write(f"Metal Flips: {sim.move_stats['accepted']['flip']}\n")
        f.write(f"Insertions: {sim.move_stats['accepted']['insert']}\n")
        f.write(f"Deletions: {sim.move_stats['accepted']['delete']}\n")
        f.write(f"1NN Shuffles: {sim.move_stats['accepted']['shuffle']}\n")
        f.write(f"MD Simulations: {sim.move_stats['accepted']['MD']}\n\n")
        
        f.write(f'Rejected Moves\n')
        f.write(f"Metal Swaps: {sim.move_stats['rejected']['swap']}\n")
        f.write(f"Interstitial Swaps: {sim.move_stats['rejected']['swap_ints']}\n")
        f.write(f"New Hosts: {sim.move_stats['rejected']['new_host']}\n")
        f.write(f"Metal Flips: {sim.move_stats['rejected']['flip']}\n")
        f.write(f"Insertions: {sim.move_stats['rejected']['insert']}\n")
        f.write(f"Deletions: {sim.move_stats['rejected']['delete']}\n")
        f.write(f"1NN Shuffles: {sim.move_stats['rejected']['shuffle']}\n")
        f.write(f"MD Simulations: {sim.move_stats['rejected']['MD']}\n")
        
        f.write(f'\nSteps for MD: {sim.md_step_counter}')

def parse_mc_statistics(file_path='data/MonteCarloStatistics'):
    """
    Re-implementation of the original parsing logic to resume simulations.
    """
    accepted_counts = {}
    rejected_counts = {}
    save = accepted = rejected = steps_completed = 0
    current_section = None

    move_key_map = {
        'Metal Swaps': 'swap',
        'Interstitial Swaps': 'swap_ints',
        'New Hosts': 'new_host',
        'Metal Flips': 'flip',
        'Insertions': 'insert',
        'Deletions': 'delete',
        '1NN Shuffles': 'shuffle',
        'MD Simulations': 'MD'
    }

    if not os.path.exists(file_path):
        return None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line: continue

            if line == 'Accepted Moves':
                current_section = 'accepted'
                continue
            elif line == 'Rejected Moves':
                current_section = 'rejected'
                continue

            if ':' in line:
                key, value = line.split(':')
                key, value = key.strip(), float(value.strip())

                if key == 'Steps Completed':
                    steps_completed = int(value)
                elif key == 'Steps for MD':
                    save = int(value)
                elif key in move_key_map:
                    move = move_key_map[key]
                    if current_section == 'accepted':
                        accepted_counts[move] = int(value)
                    elif current_section == 'rejected':
                        rejected_counts[move] = int(value)

    return {
        'steps_completed': steps_completed,
        'md_step_counter': save,
        'accepted_counts': accepted_counts,
        'rejected_counts': rejected_counts
    }