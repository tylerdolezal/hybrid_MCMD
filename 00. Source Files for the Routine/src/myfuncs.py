from ase.neighborlist import natural_cutoffs
import src.lammps_functions as fun
from ase.io import read, write
from ase.build import bulk
import pandas as pd
import numpy as np
import random
import copy
import os

def update_phase_field_dataset(PFM_data, old_system, system, dE, swap_pairs, species_counts, move_type):

    if move_type == 'shuffle':
        with open('shuffle_pairs') as file:
            line = file.readline().strip()  # Read the first line and strip any extra whitespace/newline characters
            parts = line.split(',')  # Split the line at the comma to separate the indices
            swap_pairs = [tuple(map(int, parts))]

    if move_type == 'cluster_hop':
        with open('cluster_hop_pairs') as file:
            swap_pairs = []
            for line in file:
                line = line.strip()  # Strip any extra whitespace/newline characters
                parts = line.split(',')  # Split the line at the comma to separate the indices
                swap_pairs.append(tuple(map(int, parts)))



    for pair in swap_pairs:

        atom_index1, atom_index2 = pair[0], pair[1]

        if system[atom_index2].symbol not in ['B', 'C', 'O', 'N', 'H']:
            # Get initial and final positions for both atoms
            atom_type1 = system[atom_index1].symbol
            species_counts[atom_type1] += 1
            r0_1 = old_system[atom_index1].position  # Position before the swap
            rf_1 = system[atom_index1].position  # Position after the swap, i.e., the position of the second atom

            # Log the swap for the first atom
            new_row1 = pd.DataFrame({'Atom Type': [atom_type1], 'r_0': [r0_1], 'r_f': [rf_1], 'Times Moved': [species_counts[atom_type1]], 'dE': [dE]})
            PFM_data = pd.concat([PFM_data, new_row1], ignore_index=True)

        atom_type2 = system[atom_index2].symbol
        species_counts[atom_type2] += 1
        r0_2 = old_system[atom_index2].position  # Position before the swap
        rf_2 = system[atom_index2].position  # Position after the swap, i.e., the position of the first atom

        # Log the swap for the second atom
        new_row2 = pd.DataFrame({'Atom Type': [atom_type2], 'r_0': [r0_2], 'r_f': [rf_2], 'Times Moved': [species_counts[atom_type2]], 'dE': [dE]})
        PFM_data = pd.concat([PFM_data, new_row2], ignore_index=True)

    PFM_data.to_csv('data/phase_field_model_data.csv', index=False)
    return PFM_data, species_counts


def initialize_mc_step_counter(atoms):
    # Get unique chemical symbols, sorted alphabetically
    unique_symbols = sorted(set(atoms.get_chemical_symbols()))

    # Initialize dictionary to track MC steps for each species
    mc_steps_per_species = {symbol: 0 for symbol in unique_symbols}

    return mc_steps_per_species


def snapshots(mc_step, snapshot_every):
    if mc_step % snapshot_every == 0:
        os.system(f'cp data/MonteCarloStatistics data/epochs/mc_stats_{mc_step}')
        os.system(f'cp data/species_counter.json data/epochs/species_counter_{mc_step}.json')


def place_near_host(atoms, host_index, bc_index, cutoff=2.25):
    host_position = atoms[host_index].position
    bc_position = atoms[bc_index].position
    min_distance = 2.0 # angstroms from itself in the current interstitial site
    distance = 1.0  # angstroms from the host
    sqrt_distance = 1.2*np.sqrt(2)

    # Displacement vectors to place B/C atom diagonally from the host
    displacement_vectors = [
        np.array([sqrt_distance, 0.0, 0.0]),
        np.array([-sqrt_distance, 0.0, 0.0]),
        np.array([0.0, sqrt_distance, 0.0]),
        np.array([0.0, -sqrt_distance, 0.0]),
        np.array([distance, distance, 0.0]),
        np.array([distance, -distance, 0.0]),
        np.array([-distance, distance, 0.0]),
        np.array([-distance, -distance, 0.0])]

    '''
    # Try to place the B/C atom at a nearby position
    for displacement in displacement_vectors:
        new_position = host_position + displacement

        # Ensure the new position is not the same as the current dopant position
        if not np.allclose(new_position, bc_position, atol=min_distance):  # Explicit check for dopant atom
            # Also ensure the new position is not too close to any other atoms
            if not any(np.allclose(new_position, atom.position, atol=min_distance-1.0) for atom in atoms):
                return new_position  # Return the new position if valid
    '''
    # If all initial positions are occupied, try additional spots around nearest neighbors

    # Set up neighbor list to find neighbors within the cutoff
    cutoff = natural_cutoffs(atoms)
    neighbor_list = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
    neighbor_list.update(atoms)

    # Find indices of neighbors
    indices, offsets = neighbor_list.get_neighbors(host_index)

    # Collect the types of these neighbors that are metal types
    metal_neighbors = [idx for idx in indices if atoms[idx].symbol not in ['B', 'C', 'O', 'N', 'H']]
    for attempts in range(100):
        host_index = random.choice(metal_neighbors)

        host_position = atoms[host_index].position
        # Try to place the B/C atom at a nearby position
        for displacement in displacement_vectors:
            new_position = host_position + displacement

            # Ensure the new position is not the same as the current dopant position
            if not np.allclose(new_position, bc_position, atol=min_distance):  # Explicit check for dopant atom
                # Also ensure the new position is not too close to any other atoms
                if not any(np.allclose(new_position, atom.position, atol=min_distance-1.0) for atom in atoms):
                    return new_position  # Return the new position if valid

    # not possible to place nearby
    return None

from ase import neighborlist
'''
def perform_cluster_move(system, cutoff=2.5, max_displacement=0.5):
    # Identify all B atoms
    b_indices = [atom.index for atom in system if atom.symbol in ['B', 'C', 'O', 'N', 'H']]

    # Select a random B atom
    b_index = random.choice(b_indices)
    b_position = system[b_index].position

    # Find neighbors within the cutoff distance
    neighbor_indices = []
    nl = neighborlist.NeighborList([cutoff / 2.0] * len(system), self_interaction=False, bothways=True)
    nl.update(system)

    indices, offsets = nl.get_neighbors(b_index)

    # Include the B atom itself in the cluster
    cluster_indices = [b_index] + indices

    # Apply a random displacement to the entire cluster
    displacement = (np.random.random(3) - 0.5) * max_displacement
    for index in cluster_indices:
        system[index].position += displacement

    return system
'''
from scipy.spatial.transform import Rotation as R
def perform_cluster_move(system, cutoff=2.5, max_displacement=0.5):
    # Identify all B atoms
    b_indices = [atom.index for atom in system if atom.symbol in ['B', 'C', 'O', 'N', 'H']]

    # Select a random B atom
    b_index = random.choice(b_indices)

    # Find neighbors within the cutoff distance
    cutoff = natural_cutoffs(system)

    nl = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
    nl.update(system)

    indices, offsets = nl.get_neighbors(b_index)

    # Include the B atom itself in the cluster
    cluster_indices = [b_index] + indices

    # Randomly generate a rotation axis (normalized)
    rotation_axis = np.random.normal(size=3)
    rotation_axis /= np.linalg.norm(rotation_axis)

    # Randomly select a rotation angle between 10 and 20 degrees
    rotation_angle = np.radians(random.uniform(20, 30))

    # Create rotation object
    rotation = R.from_rotvec(rotation_angle * rotation_axis)

    # Get the position of the center atom
    center_position = system.positions[b_index]

    # Apply rotation to each atom in the cluster
    for idx in cluster_indices:
        if idx != b_index:
            relative_position = system.positions[idx] - center_position
            rotated_position = rotation.apply(relative_position)
            system.positions[idx] = rotated_position + center_position

    return system

def perform_uniform_cluster_hop(system, r_cutoff=2.25, hop_distance=2.0):
    # Identify all B or C atoms
    bc_indices = [atom.index for atom in system if atom.symbol in ['B', 'C', 'O', 'N', 'H']]

    # Select a random B or C atom to move
    b_index = random.choice(bc_indices)

    # Find neighbors within the cutoff distance of the move atom (these will be metal atoms)
    b_position = system[b_index].position
    nl = neighborlist.NeighborList([r_cutoff / 2.0] * len(system), self_interaction=False, bothways=True)
    nl.update(system)

    indices, offsets = nl.get_neighbors(b_index)
    metal_indices = [
        index for index, offset in zip(indices, offsets)
        if system[index].symbol not in ['B', 'C', 'O', 'N', 'H'] and np.linalg.norm(system[index].position + np.dot(offset, system.get_cell()) - b_position) <= r_cutoff
    ]

    nonmetal_indices = [
        index for index, offset in zip(indices, offsets)
        if system[index].symbol in ['B', 'C', 'O', 'N', 'H'] and np.linalg.norm(system[index].position + np.dot(offset, system.get_cell()) - b_position) <= r_cutoff
    ]


    # Determine a random cardinal direction for the hop
    directions = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
    direction = random.choice(directions)

    # Calculate the new positions for the metal atoms
    new_positions = [system[i].position + hop_distance * direction for i in metal_indices]

    # Identify the nearest metal atoms to these new positions and swap their positions with the original cluster atoms
    nonmetal_companions = {}
    swapped_pairs = []
    displacements = []
    for old_index, new_position in zip(metal_indices, new_positions):
        distances = np.linalg.norm(system.get_positions() - new_position, axis=1)
        sorted_indices = np.argsort(distances)

        # Find the nearest metal atom to swap with
        for nearest_index in sorted_indices:
            if system[nearest_index].symbol not in ['B', 'C', 'O', 'N', 'H']:
                break

        # If a nonmetal companion is found, keep track of it
        if system[nearest_index].symbol in ['B', 'C', 'O', 'N', 'H']:
            nonmetal_companions[nearest_index] = old_index

        # Swap positions if a suitable nearest metal atom is found
        if system[nearest_index].symbol not in ['B', 'C', 'O', 'N', 'H']:
            temp_position = system[old_index].position.copy()

            displacements.append(system[nearest_index].position - system[old_index].position)

            system[old_index].position = system[nearest_index].position
            system[nearest_index].position = temp_position
            swapped_pairs.append((old_index, nearest_index))


    # Calculate the displacement for the B or C atom based on the movement of its nearest host(s)
    displacement = np.mean(displacements)

    # Move the nonmetal companions along with their metal hosts
    for nonmetal_index, host_index in nonmetal_companions.items():
        host_displacement = -displacement
        system[nonmetal_index].position += host_displacement
        swapped_pairs.append((0, nonmetal_index))

    # Move the B or C neighbors to the new cluster location
    for index in nonmetal_indices:
        system[index].position += displacement
        swapped_pairs.append((0, index))

    # Move the B or C host to the new cluster location
    system[b_index].position += displacement
    swapped_pairs.append((0, b_index))

    with open('cluster_hop_pairs', 'w') as file:
        for old_index, nearest_index in swapped_pairs:
            file.write(f"{old_index}, {nearest_index}\n")


    return system



def shuffle_neighbor_types(system, cutoff=2.25):
    """
    Selects a random B atom, finds its neighbors within a specified cutoff,
    and shuffles the atomic types among these neighbors.

    Parameters:
    - system (ASE Atoms object): The atomic configuration for simulation.
    - atom_type (str): The type of the central atom around which neighbors are identified.
    - cutoff (float): The radius within which to find neighbors, in angstroms.

    Returns:
    - system (ASE Atoms object): The modified atomic configuration after the shuffle.
    """

    b_indices = [atom.index for atom in system if atom.symbol in ['B', 'C', 'O', 'N', 'H']]

    if not b_indices:
        return system  # Return unchanged if no B atoms found

    # Select a random B atom
    b_index = random.choice(b_indices)

    # Set up neighbor list to find neighbors within the cutoff
    cutoff = natural_cutoffs(system)
    neighbor_list = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
    neighbor_list.update(system)

    # Find indices of neighbors
    indices, offsets = neighbor_list.get_neighbors(b_index)

    # Collect the types of these neighbors that are metal types
    metal_neighbors = [idx for idx in indices if system[idx].symbol not in ['B', 'C', 'O', 'N', 'H']]

    # select a nearest metal neighbor
    neighbor_index = random.choice(metal_neighbors)

    indices, offsets = neighbor_list.get_neighbors(neighbor_index)


    # Exclude the original B/C atom from the list of potential switch candidates
    indices = [idx for idx in indices if system[idx].symbol not in ['B', 'C', 'O', 'N', 'H']]

    # filter out neighbors of same type
    indices = [idx for idx in indices if system[idx].symbol != system[neighbor_index].symbol]

    if indices:
        switch_with_index = random.choice(indices)
        # Switch positions
        p1 = system.positions[neighbor_index].copy()
        p2 = system.positions[switch_with_index].copy()

        system.positions[neighbor_index], system.positions[switch_with_index] = p2, p1
        with open('shuffle_pairs', 'w') as file:
            file.write(f"{neighbor_index}, {switch_with_index}")

    # if we are in a region where no nearest neighbors of a different type exist,
    # then grab a different type metal from somewhere in the cell and try it
    else:

        # Set up neighbor list to find neighbors within a wider cutoff
        neighbor_list = neighborlist.NeighborList([4.0/2] * len(system), self_interaction=False, bothways=True)
        neighbor_list.update(system)

        indices, offsets = neighbor_list.get_neighbors(neighbor_index)

        # Exclude the original B/C atom from the list of potential switch candidates
        neighbor_indices = [idx for idx in indices if system[idx].symbol not in ['B', 'C', 'O', 'N', 'H']]

        # filter out neighbors of same type
        neighbor_indices = [idx for idx in indices if system[idx].symbol != system[neighbor_index].symbol]

        if neighbor_indices:
            switch_with_index = random.choice(neighbor_indices)
            # Switch positions
            p1 = system.positions[neighbor_index].copy()
            p2 = system.positions[switch_with_index].copy()

            system.positions[neighbor_index], system.positions[switch_with_index] = p2, p1
            with open('shuffle_pairs', 'w') as file:
                file.write(f"{neighbor_index}, {switch_with_index}")


    return system



# Relax the pre and post swapped cells
def calculate_energy_change(system, energy, swapped_pairs, move_type, run_MD, supcomp_command):

    system_copy = copy.deepcopy(system)

    if move_type == 'swap':

        for apair in swapped_pairs:
            p1 = system_copy[apair[0]].position.copy()
            p2 = system_copy[apair[1]].position.copy()

            system_copy[apair[0]].position, system_copy[apair[1]].position = p2, p1

    elif move_type == 'translate':

        for apair in swapped_pairs:

            displacement = np.random.uniform(-0.5, 0.5, size=(1, 3))

            p1 = system_copy[apair[0]].position.copy() + displacement
            p2 = system_copy[apair[1]].position.copy() + displacement

            system_copy[apair[0]].position, system_copy[apair[1]].position = p1, p2

    elif move_type == 'new_host':

        for apair in swapped_pairs:

            # introduce the B/C atom near the host
            p2 = place_near_host(system_copy, apair[0], apair[1])

            system_copy[apair[1]].position = p2

    elif move_type == 'cluster':


        system_copy = perform_cluster_move(system_copy)

    elif move_type == 'shuffle':


        system_copy = shuffle_neighbor_types(system_copy)

    elif move_type == 'cluster_hop':

        system_copy = perform_uniform_cluster_hop(system_copy)

    elif move_type == 'MD':
        system_copy, new_energy = run_md_simulation(system_copy)

    original_energy = energy[-1]

    if not run_MD:
        system_copy, new_energy = relax(system_copy, move_type, supcomp_command)

    delta_E = new_energy - original_energy

    return delta_E, new_energy

# logic for relaxing the swapped cells
def relax(system, move_type, supcomp_command):

    write("POSCAR", system, format="vasp", direct=True, sort=True)

    # create the POSCAR.data file
    fun.poscar_to_lammps()

    # execute atomic minimization without letting
    # the simcell relax from NPT sims

    if move_type == 'new_host':
        os.system(f"{supcomp_command} -in relax.in")
    else:
        os.system(f"{supcomp_command} -in static_relax.in")

    # create the CONTCAR file
    fun.lammps_to_poscar()

    system = read("CONTCAR", format='vasp')

    E = np.loadtxt("energy.txt")

    return system, E

# logic for relaxing the swapped cells
def initial_relax(system, supcomp_command):

    write("POSCAR", system, format="vasp", direct=True, sort=True)

    # create the POSCAR.data file
    fun.poscar_to_lammps()

    # execute the CG minimization simulation
    # plus the initial NPT MD
    os.system(f"{supcomp_command} -in relax-0.in")

    # create the CONTCAR file
    fun.lammps_to_poscar()

    system = read("CONTCAR", format='vasp')

    E = np.loadtxt("energy.txt")

    return system, E

# logic for executing the MD simulation
def run_md_simulation(system, supcomp_command):
    # Setting up MD with the NVT ensemble

    write("POSCAR", system, format="vasp", direct=True)

    # create the POSCAR.data file
    fun.poscar_to_lammps()

    # run NVT MD sim
    os.system(f"{supcomp_command} -in MD.in")

    # create the CONTCAR file
    fun.lammps_to_poscar()

    system = read("CONTCAR", format='vasp')

    E = np.loadtxt("energy.txt")

    return system, E

from ase import Atoms
from math import floor
def place_additives_nearby(in625_supercell, additives, O2, size):
    parent = additives[0]
    minority = additives[1]
    num_tib2_units = additives[2]

    Natoms = len(in625_supercell)

    nO2 = max(1,int(floor(0.5*num_tib2_units)))

    # Randomly select positions in the supercell to replace with TiB2
    positions = np.random.choice(range(Natoms), num_tib2_units+nO2, replace=False)

    # place our ceramics
    for pos in positions[:num_tib2_units]:
        #in625_supercell[pos].symbol = parent
        # Add B atoms close to Ti in a realistic manner
        in625_supercell += Atoms(minority, positions=[in625_supercell[pos].position + [1.25, 0, 0]])

    if O2:
        # place the O2 if we are considering O2
        for pos in positions[num_tib2_units:]:
            in625_supercell[pos].symbol = parent
            # Add B atoms close to Ti in a realistic manner
            in625_supercell += Atoms('O', positions=[in625_supercell[pos].position + [0, 1.25, 0]])



class AtomsWithVacancies(Atoms):
    def __init__(self, symbols, positions, cell, vacancies=None, **kwargs):
        super().__init__(symbols=symbols, positions=positions, cell=cell, **kwargs)
        self.vacancies = vacancies if vacancies is not None else []

    # spread vacancies out
    def introduce_vacancies(self):
        # Count the number of 'B' and 'C' atoms
        symbols = self.get_chemical_symbols()
        b_c_count = sum(1 for s in symbols if s in ['B', 'C', 'N', 'H'])

        if b_c_count > 0:
            # Calculate the number of vacancies to introduce
            num_vacancies = int(np.round(b_c_count / 2))
        else:
            # If no 'B' or 'C' atoms are present, set vacancies to 1% of total atoms
            num_vacancies = max(1, int(np.round(len(symbols) * 0.01)))

        # Determine the majority species
        unique, counts = np.unique(symbols, return_counts=True)
        majority_species = unique[np.argmax(counts)]

        # Identify indices of majority species atoms
        majority_indices = [i for i, symbol in enumerate(symbols) if symbol == majority_species]

        # Get positions of majority species atoms
        positions = self.get_positions()[majority_indices]

        # Randomly select atoms to delete to introduce vacancies
        indices_to_delete = [np.random.choice(majority_indices)]

        for _ in range(1, num_vacancies):
            distances = np.zeros((len(majority_indices), len(indices_to_delete)))

            for i, idx in enumerate(majority_indices):
                for j, del_idx in enumerate(indices_to_delete):
                    distances[i, j] = np.linalg.norm(positions[i] - self[del_idx].position)

            min_distances = distances.min(axis=1)
            next_idx = majority_indices[np.argmax(min_distances)]
            indices_to_delete.append(next_idx)

        self.vacancies.extend(indices_to_delete)
        del self[indices_to_delete]



    def get_vacancies(self):
        return self.vacancies



# Initialization
def initialize_system(composition, grain_boundary, supcomp_command, md_params, additives, O2, crystal_shape, size, vacancies):

    # set parent to most prominent species
    parent = max(composition, key=lambda x: max(composition))

    # Create a unit cell based on the parent constituents pure self
    if not grain_boundary:
        try:
            a = bulk(parent, cubic=True).cell[0,0]
            atoms = bulk('X', a = a, crystalstructure=crystal_shape, cubic=True) ## the parent lattice
        except:
            a = bulk(parent, orthorhombic=True).cell[0,0]
            atoms = bulk('X', a = a, crystalstructure=crystal_shape, orthorhombic=True) ## the parent lattice

        atoms = atoms*(size,size,size)


        # Calculate the total number of atoms
        total_atoms = len(atoms)

        # Extract symbols and percentages from the library
        symbols = list(dict(composition).keys())
        percentages = [x for x in composition.values()]

        # Calculate the exact number of each type of atom
        num_atoms = {symbol: int(round(total_atoms * percentage)) for symbol, percentage in zip(symbols, percentages)}

        # Adjust for rounding errors by adding/removing atoms from the most frequent species
        actual_total = sum(num_atoms.values())
        difference = total_atoms - actual_total
        most_frequent_symbol = symbols[np.argmax(percentages)]
        num_atoms[most_frequent_symbol] += difference


        # Create the list of atomic symbols based on the exact numbers
        atomic_symbols = []
        for symbol, count in num_atoms.items():
            atomic_symbols.extend([symbol] * count)

        # Shuffle the list to randomize the positions
        np.random.shuffle(atomic_symbols)

        # Update the symbols of the existing atoms object with atomic symbols
        atoms.set_chemical_symbols(atomic_symbols)

    else:
        # load in the pre-existing grain boundary structure
        atoms = read("POSCAR-gb", format="vasp")

    if additives:
        place_additives_nearby(atoms, additives, O2, size)

    if vacancies:
        atoms = AtomsWithVacancies(symbols=atoms.get_chemical_symbols(), positions=atoms.get_positions(), cell=atoms.get_cell())
        atoms.introduce_vacancies()


    # this way, the modfiles can be updated properly
    write("POSCAR-0", atoms, format='vasp', direct=True, sort=True)
    # update modfiles for the relaxations and md runs
    fun.update_modfiles(md_params)

    atoms, _ = initial_relax(atoms, supcomp_command)

    # save the relaxed structure as our initial structure
    write("POSCAR-0", atoms, format='vasp', direct=True, sort=True)

    # save the relaxed structure as current structure
    write("POSCAR-1", atoms, format='vasp', direct=True, sort=True)

    return(atoms)


from ase.neighborlist import NeighborList

def get_nearest_neighbors(system, atom_index, disperse=False, cutoff=2.25):
    cutoff = natural_cutoffs(system)

    if disperse:
        cutoff = [2.8 / 2.0] * len(system)

    neighbor_list = NeighborList(cutoff, self_interaction=False, bothways=True)
    neighbor_list.update(system)
    indices, offsets = neighbor_list.get_neighbors(atom_index)

    host_position = system[atom_index].position

    distances = []
    for index, offset in zip(indices, offsets):
        d = np.linalg.norm(system[index].position + np.dot(offset, system.get_cell()) - host_position)
        distances.append(d)

    # if we are attempting to disperse then do not filter out
    # the interstitial neighbors
    if disperse:
        return(indices)

    # Identify the current nearest neighbor
    nearest_neighbor_index = indices[np.argmin(distances)]

    # Collect the types of these neighbors that are metal types, omitting 'B', 'C', 'O', 'N', and the nearest neighbor
    metal_neighbors = [
        idx for idx in indices
        if system[idx].symbol not in ['B', 'C', 'O', 'N', 'H'] and idx != nearest_neighbor_index
    ]

    return metal_neighbors


def select_random_atoms(system, move_type):

    # filter out B and C atoms as we will find new 'hosts' for them rather than attempt
    # swaps and translational moves
    all_metal_indices = [i for i, atom in enumerate(system) if atom.symbol not in ['B', 'C', 'O', 'N', 'H']]
    bc_indices = [i for i, atom in enumerate(system) if atom.symbol in ['B', 'C', 'O', 'N', 'H']]

    # cannot be an odd int because these are pairs of atoms
    selection = {     'swap': 2,
                 'translate': 2,
                  'new_host': 2}


    # Perturb or swap hosts for 1 pair of atoms (2 atoms)
    num_atoms_to_swap = selection[move_type]

    if move_type != 'new_host':

        swapped = []
        while len(swapped) < num_atoms_to_swap:
            # Choose the host atom randomly from metal indices
            host_atom = random.choice(all_metal_indices)
            host_atom_type = system[host_atom].symbol  # Get the chemical symbol of the host atom
            swapped.append(host_atom)

            # Now select a nearest neighbor that is NOT of the same chemical type
            neighbors = get_nearest_neighbors(system, host_atom)
            valid_neighbors = [atom for atom in neighbors if system[atom].symbol != host_atom_type]

            if valid_neighbors:
                atom = random.choice(valid_neighbors)
                if atom not in swapped:
                    swapped.append(atom)
            else:

                # start checking 2NN shell
                neighbors = get_nearest_neighbors(system, host_atom, 4.0)
                valid_neighbors = [atom for atom in neighbors if system[atom].symbol != host_atom_type]

                if valid_neighbors:
                    atom = random.choice(valid_neighbors)
                    if atom not in swapped:
                        swapped.append(atom)


        # Pair the atoms for swapping
        swap_pairs = [(swapped[i], swapped[i+1]) for i in range(0, num_atoms_to_swap, 2)]


    elif move_type == 'new_host':

        atoms_picked = []
        attempts = 0

        choice = 0
        # do not attempt cluster growth until
        # enough time has been allowed for B/C diffusion
        if len(np.loadtxt('data/energies')) >= 500:
            choice = random.choice([0, 0, 1, 2])

            # we need to catch an error that occurs when we only have
            # a single interstitial floating around
            if len(bc_indices) < 2:
                choice = 0

        while len(atoms_picked) < num_atoms_to_swap:

            # first we try diffusing in local area
            # if no sites are available, explore the whole cell
            if attempts > 20:
                non_bc_indices = [i for i, atom in enumerate(system) if atom.symbol not in ['B', 'C', 'O', 'N', 'H']]

                # too crowded to place nearby, so flip choice back to 0
                choice = 0

            b_or_c = random.choice(bc_indices)
            non_bc_indices = get_nearest_neighbors(system, b_or_c)

            # introduce B into another B's cluster to attempt growth of phase
            if choice == 1:




                first_b_or_c = random.choice(bc_indices)

                non_bc_indices = get_nearest_neighbors(system, first_b_or_c)

                # do not grab the same B/C atom
                b_or_c = random.choice([x for x in bc_indices if x != first_b_or_c])

                metal_host = random.choice(non_bc_indices)

            # Attempt to disperse clusters
            elif choice == 2:
                # Detect if there's a cluster
                cluster_detected = False
                cluster_pair = None

                for interstitial in bc_indices:
                    neighbors = get_nearest_neighbors(system, interstitial, True)
                    neighbor_interstitials = [n for n in neighbors if system[n].symbol in ['B', 'C', 'O', 'N', 'H']]

                    if neighbor_interstitials:  # A cluster exists
                        cluster_detected = True
                        cluster_pair = (interstitial, random.choice(neighbor_interstitials))
                        break

                if cluster_detected:


                    # Select the interstitial to move
                    b_or_c = cluster_pair[0]

                    # Find distant sites outside the nearest neighbor shell
                    cluster_neighbors = get_nearest_neighbors(system, b_or_c)
                    potential_sites = [i for i in non_bc_indices if i not in cluster_neighbors]

                    if potential_sites:
                        # Select a distant site
                        metal_host = random.choice(potential_sites)
                    else:
                        # Fallback to a random non-clustered site if no distant site exists
                        metal_host = random.choice(non_bc_indices)
                else:
                    # If no clusters are detected, fallback to general diffusion
                    metal_host = random.choice(non_bc_indices)
            else:
                # Default behavior
                metal_host = random.choice(non_bc_indices)


            if metal_host not in atoms_picked and b_or_c not in atoms_picked:
                atoms_picked.append(metal_host)
                atoms_picked.append(b_or_c)

            attempts += 1

        # Pair the atoms for host swap
        swap_pairs = [(atoms_picked[i], atoms_picked[i+1]) for i in range(0, num_atoms_to_swap, 2)]

    return swap_pairs


def parse_mc_statistics(file_path='data/MonteCarloStatistics'):
    stats = {}


    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Split line into key and value
            if ':' in line:
                key, value = line.split(':')
                key = key.strip()
                value = value.strip()
                # Convert value to appropriate type
                if key in ['Steps Completed', 'Accepted Swaps', 'New Hosts Accepted', 'Translates Accepted', 'Cluster Hops Accepted', 'Cluster Shuffles Accepted', 'MD Simulations Accepted', 'Steps for MD']:
                    stats[key] = int(value)
                elif key in ['Acceptance %', 'Rejection %']:
                    stats[key] = float(value)

    move_counts = {'swap': stats['Accepted Swaps'],
                   'new_host': stats['New Hosts Accepted'],
                   'translate': stats['Translates Accepted'],
                   'cluster_hop': stats['Cluster Hops Accepted'],
                   'shuffle': stats['Cluster Shuffles Accepted'],
                   'MD': stats['MD Simulations Accepted']}

    steps_completed = stats['Steps Completed']
    acceptance_count = int((stats['Acceptance %'] / 100) * steps_completed)
    rejection_count = int((stats['Rejection %'] / 100) * steps_completed)

    return stats['Steps for MD'], acceptance_count, rejection_count, move_counts, steps_completed
