from ase.neighborlist import natural_cutoffs
import src.lammps_functions as fun
from ase.io import read, write
from ase.build import bulk
import pandas as pd
import numpy as np
import random
import copy
import os

interstitials = ["B", "C", "Cl", "O", "N", "H"]
ignore = []

common_structures = {
    # Nickel-based superalloy elements
    "Ni": "fcc",   # Base element for Ni superalloys
    "Co": "hcp",   # Strengthens γ' phase (can be fcc at high temp)
    "Cr": "bcc",   # Oxidation resistance, stabilizes γ phase
    "Fe": "bcc",   # Present in some superalloys & stainless steels
    "Mo": "bcc",   # Strengthens γ phase, creep resistance
    "Nb": "bcc",   # Forms NbC/Nb3C phases for carbide strengthening
    "Ti": "hcp",   # Strengthens γ' phase (TiAl)
    "Al": "fcc",   # Key for γ' phase formation (Ni3Al)
    "W": "bcc",    # Solid solution strengthening, creep resistance
    "Ta": "bcc",   # Similar role as W, improves oxidation resistance
    "Re": "hcp",   # Used in advanced single-crystal superalloys
    "Hf": "hcp",   # Strengthens grain boundaries in SX alloys

    # Common stainless steel elements
    "Mn": "bcc",   # Affects work-hardening and toughness
    "Si": "diamond",  # Improves oxidation resistance, aids deoxidation
    "V": "bcc",    # Strengthens carbide phases
    "Cu": "fcc",   # Sometimes added for corrosion resistance
}
# Global freeze_threshold (default to no freezing)
freeze_threshold = 0.0

def set_global_threshold(surface):
    """Sets the global freeze_threshold based on the surface configuration."""
    global freeze_threshold  # Declare global before modifying it

    if surface and isinstance(surface, list) and len(surface) > 2:
        freeze_threshold = float(surface[2])  # Extract height value safely
    else:
        freeze_threshold = 0.0  # Default: No freezing
    

    print(f"Global freeze_threshold set to: {freeze_threshold}")
    return freeze_threshold

def set_interstitials(additives):
    global interstitials
    interstitials = list(set(interstitials) | set(additives))  # Merge and remove duplicates

from ase.neighborlist import NeighborList
class MaskedNeighborList(NeighborList):
    def __init__(self, cutoff, system, **kwargs):
        """
        A universal NeighborList wrapper that ignores frozen atoms in lower layers.
        
        Parameters:
            cutoff (float): Cutoff radius for nearest neighbor search.
            system (ase.Atoms): The simulation cell.
            kwargs: Additional arguments passed to NeighborList.
        """
        global freeze_threshold  # Use the global setting

        self.system = system
        self.freeze_threshold = freeze_threshold  # Set once and never change
        
        # Identify frozen atoms (z ≤ threshold)
        self.frozen_mask = system.positions[:, 2] <= self.freeze_threshold
        
        # Initialize base NeighborList
        super().__init__(cutoff, self_interaction=False, bothways=True, **kwargs)
        self.update(system)  # Update with the initial system

    def get_neighbors(self, atom_index):
        """Overrides get_neighbors() to exclude frozen atoms."""
        if self.frozen_mask[atom_index]:
            return [], []  # No neighbors for frozen atoms

        indices, offsets = super().get_neighbors(atom_index)
        
        # Filter out frozen neighbors
        valid_neighbors = [idx for idx in indices if not self.frozen_mask[idx]]
        valid_offsets = [offsets[i] for i in range(len(indices)) if not self.frozen_mask[indices[i]]]

        return valid_neighbors, valid_offsets


def update_phase_field_dataset(PFM_data, old_system, system, dE, swap_pairs, species_counts, move_type):

    if move_type == 'shuffle':
        with open('shuffle_pairs') as file:
            line = file.readline().strip()  # Read the first line and strip any extra whitespace/newline characters
            parts = line.split(',')  # Split the line at the comma to separate the indices
            swap_pairs = [tuple(map(int, parts))]

    if move_type == 'flip':
        with open('flip_pairs') as file:
            line = file.readline().strip()  # Read the first line and strip any extra whitespace/newline characters
            parts = line.split(',')  # Split the line at the comma to separate the indices
            swap_pairs = [tuple(map(int, parts))]

    

    for pair in swap_pairs:

        atom_index1, atom_index2 = pair[0], pair[1]

        if system[atom_index2].symbol not in interstitials and system[atom_index1].symbol not in ignore:
            # Get initial and final positions for both atoms
            atom_type1 = system[atom_index1].symbol
            species_counts.setdefault(atom_type1, 0)
            species_counts[atom_type1] = species_counts.get(atom_type1, 0) + 1
            r0_1 = old_system[atom_index1].position  # Position before the swap
            rf_1 = system[atom_index1].position  # Position after the swap, i.e., the position of the second atom

            # Log the swap for the first atom
            new_row1 = pd.DataFrame({'Atom Type': [atom_type1], 'r_0': [r0_1], 'r_f': [rf_1], 'Times Moved': [species_counts[atom_type1]], 'dE': [dE]})
            PFM_data = pd.concat([PFM_data, new_row1], ignore_index=True)

        atom_type2 = system[atom_index2].symbol
        species_counts.setdefault(atom_type2, 0)
        species_counts[atom_type2] = species_counts.get(atom_type2, 0) + 1
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
        os.system(f'cp POSCAR-1 structures/POSCAR-{mc_step}')
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
    neighbor_list = MaskedNeighborList(cutoff, system=atoms)
    neighbor_list.update(atoms)

    # Find indices of neighbors
    indices, offsets = neighbor_list.get_neighbors(host_index)

    # Collect the types of these neighbors that are metal types
    metal_neighbors = [idx for idx in indices if atoms[idx].symbol not in interstitials+ignore and (freeze_threshold <= 0.0 or atoms[idx].position[2] > freeze_threshold)]
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

def shuffle_neighbor_types(system, cutoff=2.25, local=False):
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

    b_indices = [atom.index for atom in system if atom.symbol in interstitials]
    metal_indices = [atom.index for atom in system if atom.symbol  and (freeze_threshold <= 0.0 or system[atom].position[2] > freeze_threshold)]

    if not b_indices:
        return system  # Return unchanged if no B atoms found

    # Select a random B atom
    b_index = random.choice(b_indices)

    # find metal neighbors nearest to the interstitial
    metal_neighbors = get_nearest_neighbors(system, b_index, cutoff=2.25) 

    # select a nearest metal neighbor
    neighbor_index = random.choice(metal_neighbors)

    if local:
        # Get the indices of the nearest neighbors of the selected neighbor
        indices = get_nearest_neighbors(system, neighbor_index, cutoff=2.25)

        # filter out neighbors of same type; ignore self-swaps
        indices = [idx for idx in indices if system[idx].symbol != system[neighbor_index].symbol]
    
    else:
        indices = [idx for idx in metal_indices if system[idx].symbol != system[neighbor_index].symbol and (freeze_threshold <= 0.0 or system[idx].position[2] > freeze_threshold)]

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
        neighbor_list = MaskedNeighborList([5.0/2] * len(system), system=system)
        neighbor_list.update(system)

        indices, offsets = neighbor_list.get_neighbors(neighbor_index)

        # Exclude the original B/C atom from the list of potential switch candidates
        neighbor_indices = [idx for idx in indices if system[idx].symbol not in interstitials+ignore and (freeze_threshold <= 0.0 or system[idx].position[2] > freeze_threshold)]

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
        
        else:
            np.savetxt("Failed to identify shuffle pairs!", [])
            return system

    return system


from collections import Counter

def generate_chemical_potentials(choices, supcomp_command):
    dir = "src/chemical_potentials/"
    if not os.path.exists(dir):
        os.makedirs(dir)
    
    chemical_potentials = {}
    for metal in choices:
        filename = f"{dir}mu_{metal}.txt"
        if not os.path.exists(filename):
            try:
                atoms = bulk(metal, crystalstructure=common_structures[metal], cubic=True).repeat((2,2,2))
            except:
                atoms = bulk(metal, crystalstructure=common_structures[metal], orthorhombic=True).repeat((2,2,2))
            
            _, mu = relax(atoms, 'flip', supcomp_command)
            np.savetxt(filename, [mu/len(atoms)])

        chemical_potentials[metal] = np.loadtxt(filename)
    
    
    # Find the lowest-energy binary compounds
    main_system = read("POSCAR")
    
    # Count occurrences of each metal in main_system
    metal_counts = Counter(atom.symbol for atom in main_system if atom.symbol not in interstitials+ignore)
    # Identify the majority metal (most frequent element)
    majority_metal = max(metal_counts, key=metal_counts.get)

    binary_potentials = {}
    for metal2 in choices:
        if metal2 == majority_metal:
            continue  # Skip cases where metal2 is the majority metal

        # Ensure sorting order for consistency
        metal1, metal2 = sorted([majority_metal, metal2])
        binary_filename = f"src/chemical_potentials/mu_{metal1}_{metal2}.txt"
        poscar_path = f"src/chemical_potentials/poscars/{metal1}{metal2}_POSCAR"

        # Check if the binary potential is already computed and saved
        if os.path.exists(binary_filename):
            binary_potentials[(metal1, metal2)] = np.loadtxt(binary_filename)
            continue  # Skip recomputation

        # If POSCAR file exists, read and relax it to get binary energy
        if os.path.exists(poscar_path):
            try:
                atoms = read(poscar_path)  # Read the POSCAR file
                _, mu = relax(atoms, 'flip', supcomp_command)  # Run relaxation function
                E_binary = mu  # Extract binary energy
            except Exception as e:
                print(f"⚠ Error processing {poscar_path}: {e}")
                continue  # Skip this pair if relaxation fails
        else:
            print(f"⚠ Missing POSCAR file: {poscar_path}. Skipping {metal1}-{metal2}")
            continue  # Skip if the POSCAR file does not exist

        # Get pure metal energies
        E_A = chemical_potentials.get(metal1, 0.0)  # Default to 0 if missing
        E_B = chemical_potentials.get(metal2, 0.0)  # Default to 0 if missing

        # Get stoichiometry from the binary compound
        num_A = sum(1 for atom in atoms if atom.symbol == metal1)
        num_B = sum(1 for atom in atoms if atom.symbol == metal2)

        if num_A == 0 or num_B == 0:
            print(f"⚠ Invalid stoichiometry detected for {metal1}-{metal2}. Skipping.")
            continue  # Skip if no valid stoichiometry

        # Compute chemical potentials
        mu_A = (E_binary - (num_B * E_B)) / num_A
        mu_B = (E_binary - (num_A * E_A)) / num_B

        # Save only the non-majority metal's chemical potential
        non_majority_mu = mu_B if metal2 != majority_metal else mu_A
        np.savetxt(binary_filename, [non_majority_mu])  # Save computed value

        # Store in dictionary
        binary_potentials[(metal1, metal2)] = non_majority_mu

    binary_potentials[(majority_metal, majority_metal)] = chemical_potentials[majority_metal]
    return majority_metal, chemical_potentials, binary_potentials

def flip_atoms(system, metal_choices, supcomp_command):
    choices = ["Cr", "Fe", "Mo", "Nb", "Ni", "Ti"]
    if metal_choices:
        choices = metal_choices

    # Generate chemical potentials and identify the majority metal
    majority_metal, chemical_potentials, binary_potentials = generate_chemical_potentials(choices, supcomp_command)

    # Grab all metal indices (ignoring B, C, H, N, O)
    indices = [atom.index for atom in system if atom.symbol not in interstitials+ignore and (freeze_threshold <= 0.0 or system[atom].position[2] > freeze_threshold)]
    if not indices:
        raise ValueError("No valid metal atoms found in the system.")

    # Grab a random metal index
    random_index = random.choice(indices)
    original_symbol = system[random_index].symbol

    # Don't choose yourself
    true_choices = [x for x in choices if x != original_symbol]
    if not true_choices:
        raise ValueError(f"No valid swap choices available for {original_symbol}.")

    flip_symbol = random.choice(true_choices)

    # Ensure sorted keys for binary lookup
    ref1, ref2 = sorted([majority_metal, original_symbol])
    ref3, ref4 = sorted([majority_metal, flip_symbol])

    # Retrieve chemical potentials with fallbacks
    mu_A = binary_potentials[(ref1, ref2)]
    mu_B = binary_potentials[(ref3, ref4)]

    # If swapping to/from majority metal, enforce pure values
    if flip_symbol == majority_metal:
        mu_B = chemical_potentials[majority_metal]
    if original_symbol == majority_metal:
        mu_A = chemical_potentials[majority_metal]

    # Assign correctly: mu_i corresponds to the original atom, mu_j to the new atom
    mu_i = mu_A if original_symbol == ref2 else mu_B
    mu_j = mu_B if flip_symbol == ref4 else mu_A

    # Flip the atom type in the system
    system[random_index].symbol = flip_symbol

    # Compute chemical potential difference
    delta_mu = float(mu_j - mu_i)  # Ensure delta_mu is always a float

    # Log the swap in a file
    with open('flip_pairs', 'w') as file:
        file.write(f"{random_index}, {original_symbol} -> {flip_symbol}\n")

    return system, delta_mu



def add_new_species(atoms, counter):

    unique_species = set([atom.symbol for atom in atoms])
    for species in unique_species:
        counter.setdefault(species, 0)
    
    return(counter)

# Relax the pre and post swapped cells
def calculate_energy_change(system, energy, swapped_pairs, move_type, run_MD, supcomp_command, metal_choices=None):

    system_copy = copy.deepcopy(system)
    delta_mu = 0.0 # for when we are not flipping

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

    elif move_type == 'shuffle':
        system_copy = shuffle_neighbor_types(system_copy)

    elif move_type == 'MD':
        system_copy, new_energy = run_md_simulation(system_copy)
    
    elif move_type == 'flip':
         system_copy, delta_mu = flip_atoms(system_copy, metal_choices, supcomp_command)

    original_energy = energy[-1]
    if not run_MD:
        system_copy, new_energy = relax(system_copy, move_type, supcomp_command)

    delta_E = (new_energy - original_energy) + delta_mu

    return delta_E, new_energy

# logic for relaxing the swapped cells
def relax(system, move_type, supcomp_command):

    write("POSCAR", system, format="vasp", direct=True, sort=False)
    
    # account for new species types
    fun.update_modfiles()
    # create the POSCAR.data file
    fun.poscar_to_lammps()

    # execute atomic minimization without letting
    # the simcell relax from NPT sims

    if move_type in ['new_host', 'flip']:
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
def place_additives_nearby(in625_supercell, additives):
    global freeze_threshold  # Use the global parameter
    for additive in additives:
        parent = additive[0]
        minority = additive[1]
        num_tib2_units = additive[2]

        Natoms = len(in625_supercell)
        # Randomly select positions in the supercell to replace with TiB2
        # Get valid (non-frozen) atom indices
        non_frozen_indices = [i for i in range(Natoms) if in625_supercell[i].position[2] > freeze_threshold]
        # Randomly select valid positions
        positions = np.random.choice(non_frozen_indices, num_tib2_units, replace=False)

        # place our ceramics
        for pos in positions[:num_tib2_units]:
            # transmute our partner into the metallic host
            if parent.lower() != 'none':
                in625_supercell[pos].symbol = parent
            # Add B atoms close to Ti in a realistic manner (account for possible digit, i.e., B2)
            try:
                in625_supercell += Atoms(minority, positions=[in625_supercell[pos].position + np.array([1.25, 0, 0])])
            except Exception as e:
                placement = [np.array([1.25, 0, 0]), np.array([0, 1.25, 0])]
                for array in placement:
                    in625_supercell += Atoms(minority[0], positions=[
                    in625_supercell[pos].position + array])


# Define Van der Waals diameters for common elements (approximate values in Å)
VAN_DER_WAALS_DIAMETERS = {
    "H": 2.4,   # Hydrogen
    "He": 2.8,  # Helium
    "Li": 3.0,  # Lithium
    "Be": 3.0,  # Beryllium
    "B": 3.9,   # Boron
    "C": 3.4,   # Carbon
    "N": 3.1,   # Nitrogen
    "O": 3.0,   # Oxygen
    "F": 2.8,   # Fluorine
    "Ne": 3.0,  # Neon
    "Na": 3.7,  # Sodium
    "Mg": 3.6,  # Magnesium
    "Al": 4.2,  # Aluminum
    "Si": 4.2,  # Silicon
    "P": 4.3,   # Phosphorus
    "S": 3.6,   # Sulfur
    "Cl": 3.5,  # Chlorine
    "Ar": 3.8,  # Argon
    "K": 4.7,   # Potassium
    "Ca": 4.5,  # Calcium
    "Sc": 4.4,  # Scandium
    "Ti": 4.2,  # Titanium
    "V": 4.1,   # Vanadium
    "Cr": 4.0,  # Chromium
    "Mn": 4.0,  # Manganese
    "Fe": 4.0,  # Iron
    "Co": 4.0,  # Cobalt
    "Ni": 4.0,  # Nickel
    "Cu": 4.2,  # Copper
    "Zn": 4.1,  # Zinc
    "Ga": 4.2,  # Gallium
    "Ge": 4.3,  # Germanium
    "As": 4.3,  # Arsenic
    "Se": 4.4,  # Selenium
    "Br": 3.7,  # Bromine
    "Kr": 4.0,  # Krypton
    "Rb": 4.8,  # Rubidium
    "Sr": 5.0,  # Strontium
    "Y": 4.8,   # Yttrium
    "Zr": 4.6,  # Zirconium
    "Nb": 4.6,  # Niobium
    "Mo": 4.5,  # Molybdenum
    "Tc": 4.4,  # Technetium
    "Ru": 4.3,  # Ruthenium
    "Rh": 4.3,  # Rhodium
    "Pd": 4.3,  # Palladium
    "Ag": 4.4,  # Silver
    "Cd": 4.6,  # Cadmium
    "In": 4.9,  # Indium
    "Sn": 5.0,  # Tin
    "Sb": 5.0,  # Antimony
    "Te": 5.1,  # Tellurium
    "I": 4.9,   # Iodine
    "Xe": 4.4,  # Xenon
    "Cs": 5.2,  # Cesium
    "Ba": 5.1,  # Barium
    "La": 5.3,  # Lanthanum
    "Ce": 5.2,  # Cerium
    "Pr": 5.2,  # Praseodymium
    "Nd": 5.2,  # Neodymium
    "Pm": 5.2,  # Promethium
    "Sm": 5.2,  # Samarium
    "Eu": 5.1,  # Europium
    "Gd": 5.2,  # Gadolinium
    "Tb": 5.2,  # Terbium
    "Dy": 5.2,  # Dysprosium
    "Ho": 5.2,  # Holmium
    "Er": 5.2,  # Erbium
    "Tm": 5.2,  # Thulium
    "Yb": 5.1,  # Ytterbium
    "Lu": 5.2,  # Lutetium
    "Hf": 4.7,  # Hafnium
    "Ta": 4.6,  # Tantalum
    "W": 4.6,   # Tungsten
    "Re": 4.6,  # Rhenium
    "Os": 4.5,  # Osmium
    "Ir": 4.4,  # Iridium
    "Pt": 4.5,  # Platinum
    "Au": 4.5,  # Gold
    "Hg": 4.7,  # Mercury
    "Tl": 5.0,  # Thallium
    "Pb": 5.2,  # Lead
    "Bi": 5.3,  # Bismuth
    "Po": 5.4,  # Polonium
    "At": 5.5,  # Astatine
    "Rn": 4.8,  # Radon
    "Fr": 5.6,  # Francium
    "Ra": 5.5,  # Radium
    "Ac": 5.5,  # Actinium
    "Th": 5.4,  # Thorium
    "Pa": 5.3,  # Protactinium
    "U": 5.2,   # Uranium
    "Np": 5.2,  # Neptunium
    "Pu": 5.2,  # Plutonium
    "Am": 5.2,  # Americium
    "Cm": 5.2,  # Curium
    "Bk": 5.1,  # Berkelium
    "Cf": 5.1,  # Californium
    "Es": 5.1,  # Einsteinium
    "Fm": 5.1,  # Fermium
    "Md": 5.1,  # Mendelevium
    "No": 5.1,  # Nobelium
    "Lr": 5.1,  # Lawrencium
}

from ase import Atoms
import warnings
def add_sheet(surface, adsorbate, lattice):
    """
    Places a sheet of an adsorbate (e.g., oxygen, chlorine) atop the uppermost layer of an ASE atoms object.
    
    Parameters:
        surface (ase.Atoms): The surface to which the adsorbate sheet will be added.
        adsorbate (str): Atomic symbol of the adsorbate (e.g., "O", "Cl").
        lattice (str): Type of lattice ("square" or "hexagonal").
        
    Returns:
        ase.Atoms: A new ASE Atoms object with the adsorbate sheet added.
    """
    vacuum_extension = 50.0
    # Get the existing atomic positions
    positions = surface.positions
    cell = surface.cell

    # Find the highest z-value (top surface)
    max_z = np.max(positions[:, 2])

    # Get the total cell height (z-dimension)
    cell_z = cell[2, 2]

    # Check if vacuum is sufficient
    available_vacuum = cell_z - max_z
    if available_vacuum < 10.0:
        vacuum_extension = vacuum_extension - available_vacuum
        warnings.warn(
            f"Insufficient vacuum above the surface (only {available_vacuum:.2f} Å). "
            f"Extending the cell by {vacuum_extension:.2f} Å."
        )
        
        cell[2, 2] += vacuum_extension
        surface.set_cell(cell)

    # Get x and y bounds from the cell with a slight margin to avoid edge effects
    x_min, x_max = np.min(positions[:, 0]) + (0.1 * cell[0][0]), np.max(positions[:, 0]) - (0.1 * cell[0][0])
    y_min, y_max = np.min(positions[:, 1]) + (0.1 * cell[1][1]), np.max(positions[:, 1]) - (0.1 * cell[1][1])

    # Determine spacing using the full van der Waals diameter
    spacing = VAN_DER_WAALS_DIAMETERS.get(adsorbate, 3.0)  # Default to 3.0 Å if unknown

    # Identify the maximum vdW diameter among the species present in the surface
    surface_species = set(surface.get_chemical_symbols())
    max_vdW_surface = max([VAN_DER_WAALS_DIAMETERS.get(species, 4.0) for species in surface_species])

    # Compute height offset using the half-sum rule
    height_offset = 0.5 * (VAN_DER_WAALS_DIAMETERS.get(adsorbate, 3.0) + max_vdW_surface)

    # Generate a 2D grid for adsorbate atoms
    x_coords = np.arange(x_min, x_max, spacing)
    y_coords = np.arange(y_min, y_max, spacing)
    
    adsorbate_positions = []
    
    if lattice == "hexagonal":
        for i, x in enumerate(x_coords):
            for j, y in enumerate(y_coords):
                shift = (spacing / 2) if i % 2 == 1 else 0  # Offset every other row
                adsorbate_positions.append([x, y + shift, max_z + height_offset])
    else:  # Default is square grid
        for x in x_coords:
            for y in y_coords:
                adsorbate_positions.append([x, y, max_z + height_offset])

    # Create an Atoms object for the adsorbate sheet
    adsorbate_layer = Atoms(symbols=[adsorbate] * len(adsorbate_positions), positions=adsorbate_positions, cell=cell, pbc=True)

    # Combine with original surface
    combined = surface + adsorbate_layer

    return combined


class AtomsWithVacancies(Atoms):
    def __init__(self, symbols, positions, cell, vacancies=None, **kwargs):
        super().__init__(symbols=symbols, positions=positions, cell=cell, **kwargs)
        self.vacancies = vacancies if vacancies is not None else []

    # spread vacancies out
    def introduce_vacancies(self):
        # Count the number of 'B' and 'C' atoms
        symbols = self.get_chemical_symbols()
        b_c_count = sum(1 for s in symbols if s in interstitials)

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
def initialize_system(composition, grain_boundary, supcomp_command, md_params, additives, crystal_shape, size, vacancies, randomize, surface):

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
        if randomize:
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

    if surface:
        atoms = add_sheet(atoms, surface[0], surface[1])

    if additives:
        place_additives_nearby(atoms, additives)

    if vacancies:
        atoms = AtomsWithVacancies(symbols=atoms.get_chemical_symbols(), positions=atoms.get_positions(), cell=atoms.get_cell())
        atoms.introduce_vacancies()


    # this way, the modfiles can be updated properly
    write("POSCAR-0", atoms, format='vasp', direct=True, sort=True)
    write("POSCAR", atoms, format='vasp', direct=True, sort=True)
    # update modfiles for the relaxations and md runs
    fun.update_modfiles(md_params)

    atoms, _ = initial_relax(atoms, supcomp_command)

    # save the relaxed structure as our initial structure
    write("POSCAR-0", atoms, format='vasp', direct=True, sort=True)

    # save the relaxed structure as current structure
    write("POSCAR-1", atoms, format='vasp', direct=True, sort=True)

    return(atoms)

def get_nearest_neighbors(system, atom_index, disperse=False, cutoff=2.25, max_cutoff=5.0):
    cutoff = natural_cutoffs(system)

    if disperse:
        cutoff = [2.8 / 2.0] * len(system)

    metal_neighbors = []

    # Incrementally increase the cutoff until we find a valid metal neighbor
    while not metal_neighbors and cutoff[0] <= max_cutoff:
        # Create a new neighbor list with the current cutoff
        neighbor_list = MaskedNeighborList(cutoff, system=system)
        neighbor_list.update(system)
        indices, offsets = neighbor_list.get_neighbors(atom_index)

        host_position = system[atom_index].position

        distances = []
        for index, offset in zip(indices, offsets):
            d = np.linalg.norm(system[index].position + np.dot(offset, system.get_cell()) - host_position)
            distances.append(d)

        # If dispersing, return all found neighbors immediately
        if disperse:
            return indices

        # Identify the current nearest neighbor
        if indices:
            nearest_neighbor_index = indices[np.argmin(distances)]
        else:
            nearest_neighbor_index = None

        # Collect only metal neighbors, ignoring interstitials & frozen atoms
        metal_neighbors = [
            idx for idx in indices
            if system[idx].symbol not in interstitials + ignore
            and idx != nearest_neighbor_index
            and (freeze_threshold <= 0.0 or system[idx].position[2] > freeze_threshold)
        ]

        # If no metal neighbors found, increase cutoff and try again
        if not metal_neighbors:
            cutoff = [min(c + 0.25, max_cutoff) for c in cutoff]

    return metal_neighbors


def select_random_atoms(system, move_type, local=False):
    global freeze_threshold  # Use the global parameter
    # filter out B and C atoms as we will find new 'hosts' for them rather than attempt
    # swaps and translational moves
    all_metal_indices = [i for i, atom in enumerate(system) if atom.symbol not in interstitials+ignore and (freeze_threshold <= 0.0 or atom.position[2] > freeze_threshold)]
    bc_indices = [i for i, atom in enumerate(system) if atom.symbol in interstitials]

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

            if local:
                # Now select a nearest neighbor that is NOT of the same chemical type
                neighbors = get_nearest_neighbors(system, host_atom)
                valid_neighbors = [atom for atom in neighbors if system[atom].symbol != host_atom_type and (freeze_threshold <= 0.0 or system[atom].position[2] > freeze_threshold)]
            

                if valid_neighbors:
                    atom = random.choice(valid_neighbors)
                    if atom not in swapped:
                        swapped.append(atom)
                else:

                    # start checking 2NN shell
                    neighbors = get_nearest_neighbors(system, host_atom, 4.0)
                    valid_neighbors = [atom for atom in neighbors if system[atom].symbol != host_atom_type and (freeze_threshold <= 0.0 or system[atom].position[2] > freeze_threshold)]

                    if valid_neighbors:
                        atom = random.choice(valid_neighbors)
                        if atom not in swapped:
                            swapped.append(atom)
            
            else:
                valid_neighbors = [atom for atom in all_metal_indices if system[atom].symbol != host_atom_type and (freeze_threshold <= 0.0 or system[atom].position[2] > freeze_threshold)]

                if valid_neighbors:
                    atom = random.choice(valid_neighbors)
                    if atom not in swapped:
                        swapped.append(atom)
                
                else:
                    np.savetxt("Failed to find valid neighbors!", [])
                    return system


        # Pair the atoms for swapping
        swap_pairs = [(swapped[i], swapped[i+1]) for i in range(0, num_atoms_to_swap, 2)]


    elif move_type == 'new_host':

        atoms_picked = []
        attempts = 0

        choice = 0
        # do not attempt cluster growth until
        # enough time has been allowed for B/C diffusion
        if len(np.loadtxt('data/energies')) >= 500:
            choice = random.choice([0, 0, 0, 1, 2])

            # we need to catch an error that occurs when we only have
            # a single interstitial floating around
            if len(bc_indices) < 2:
                choice = 0

        while len(atoms_picked) < num_atoms_to_swap:

            # first we try diffusing in local area
            # if no sites are available, explore the whole cell
            if attempts > 20:
                non_bc_indices = [i for i, atom in enumerate(system) if atom.symbol not in interstitials+ignore and (freeze_threshold <= 0.0 or atom.position[2] > freeze_threshold)]

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
                    neighbor_interstitials = [n for n in neighbors if system[n].symbol in interstitials]

                    if neighbor_interstitials:  # A cluster exists
                        cluster_detected = True
                        cluster_pair = (interstitial, random.choice(neighbor_interstitials))
                        break

                if cluster_detected:


                    # Select the interstitial to move
                    b_or_c = cluster_pair[0]

                    # Find distant sites outside the nearest neighbor shell
                    cluster_neighbors = get_nearest_neighbors(system, b_or_c)
                    potential_sites = [i for i in non_bc_indices if i not in cluster_neighbors and (freeze_threshold <= 0.0 or system[i].position[2] > freeze_threshold)]

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

            if attempts > 25:
                np.savetxt("Failed to find new host!", [])

                first = random.choice(non_bc_indices)
                atoms_picked.append(first)
                atoms_picked.append(random.choice([i for i in non_bc_indices if i not in atoms_picked and system[i].symbol != system[first].symbol]))

                # Pair the atoms for host swap
                swap_pairs = [(atoms_picked[i], atoms_picked[i+1]) for i in range(0, num_atoms_to_swap, 2)]

                return swap_pairs
                

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
                if key in ['Steps Completed', 'Accepted Swaps', 'New Hosts Accepted', 'Flips Accepted', 'Cluster Shuffles Accepted', 'MD Simulations Accepted', 'Steps for MD']:
                    stats[key] = int(value)
                elif key in ['Acceptance %', 'Rejection %']:
                    stats[key] = float(value)

    move_counts = {'swap': stats['Accepted Swaps'],
                   'new_host': stats['New Hosts Accepted'],
                   'flip': stats['Flips Accepted'],
                   'shuffle': stats['Cluster Shuffles Accepted'],
                   'MD': stats['MD Simulations Accepted']}

    steps_completed = stats['Steps Completed']
    acceptance_count = int((stats['Acceptance %'] / 100) * steps_completed)
    rejection_count = int((stats['Rejection %'] / 100) * steps_completed)

    return stats['Steps for MD'], acceptance_count, rejection_count, move_counts, steps_completed
