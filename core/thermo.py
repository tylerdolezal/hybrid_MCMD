import os
import numpy as np
import random
from ase.io import read, write
from ase import Atoms, units
from ase.optimize import FIRE, MDMin
from ase.md import Langevin
from ase.filters import UnitCellFilter
from ase.build import bulk

# --- Global Physics Constants ---
K_B = 8.617333262145e-5  # Boltzmann constant in eV/K

def initialize_system(config, calculator):
    """
    Builds the system from scratch or loads a custom cell.
    Applies initial randomization and relaxation.
    """
    comp = config['system']['composition']
    randomize = config['system']['randomize_initial']
    parent = max(comp, key=comp.get)

    if not config['system']['use_custom_cell']:
        # Build bulk lattice from scratch
        shape = config['system']['crystal_shape']
        size = config['system']['supercell_size'] # Expected as [x, y, z]

        try:
            atoms = bulk(parent, crystalstructure=shape, cubic=True)
        except:
            atoms = bulk(parent, crystalstructure=shape, orthorhombic=True)
        
        atoms = atoms * (size[0], size[1], size[2])

        if randomize:
            total_atoms = len(atoms)
            symbols = list(comp.keys())
            num_atoms = {s: int(round(total_atoms * comp[s])) for s in symbols}
            
            # Correction for rounding errors
            diff = total_atoms - sum(num_atoms.values())
            num_atoms[parent] += diff

            atomic_symbols = []
            for s, count in num_atoms.items():
                atomic_symbols.extend([s] * count)
            
            np.random.shuffle(atomic_symbols)
            atoms.set_chemical_symbols(atomic_symbols)
    else:
        # Load user-provided cell
        if not os.path.exists("POSCAR-custom"):
            raise FileNotFoundError("use_custom_cell is True but POSCAR-custom not found.")
        atoms = read("POSCAR-custom", format="vasp")

        if randomize:
            # Apply the 10-zone randomization logic for custom (e.g. GB) cells
            atoms = _apply_zonal_randomization(atoms, comp, parent)

    # Initial structure optimization
    atoms, e0 = initial_relax(atoms, calculator)
    write("POSCAR-0", atoms, format='vasp', direct=True, sort=True)
    return atoms, e0

def calculate_energy_change(system, current_energy, swap_pairs, move_type, config, calculator, mu_cache):
    """
    Evaluates probability with internalized scaling for GC and SGC ensembles.
    """
    temp = config['simulation']['temperature']
    system_copy = system.copy()
    
    delta_mu = 0.0
    mu_effective = 0.0
    prefactor = 1.0
    new_voids = None

    # --- Move Execution ---
    if move_type == 'new_host':
        idx_int, idx_host = swap_pairs[0]
        
        # Place the interstitial near the SPECIFIC host selected
        new_pos = place_near_host(system_copy, host_index=idx_host, bc_index=idx_int)
        
        if new_pos is not None:
            system_copy[idx_int].position = new_pos

        else:
            # Reject move if no geometric pocket exists around that metal atom
            return 0.0, system, current_energy, None
    
    if move_type in ['swap', 'shuffle', 'swap_ints']:
        p1, p2 = system_copy[swap_pairs[0][0]].position.copy(), system_copy[swap_pairs[0][1]].position.copy()
        system_copy[swap_pairs[0][0]].position, system_copy[swap_pairs[0][1]].position = p2, p1

    elif move_type == 'flip':
        idx = swap_pairs[0][0]
        # SGC logic: delta_mu is now passed/read from mu_cache
        delta_mu_val = mu_cache['delta_mu']

        choices = config['ensembles']['semi_grand']['metal_library']
        new_sym = random.choice([s for s in choices if s != system_copy[idx].symbol])
        system_copy[idx].symbol = new_sym
        delta_mu = -delta_mu_val if new_sym == "Ni" else delta_mu_val

    elif move_type in ['insert', 'delete']:
        # Grand Canonical Logic using cached effective mu
        system_copy, mu_add, prefactor, new_voids = insert_deletion_move_logic(
            system_copy, config['ensembles']['grand']['voids_file'], mu_cache['additives'], move_type
        )
        mu_effective = mu_add if move_type == "insert" else -mu_add

    # --- Relaxation ---
    relaxed_system, new_energy = relax_config(system_copy, move_type, calculator)

    # --- Metropolis Criterion ---
    delta_e = (new_energy - current_energy) - delta_mu - mu_effective
    probability = prefactor * np.exp(-delta_e / (K_B * temp))

    return probability, relaxed_system, new_energy, new_voids

def _apply_zonal_randomization(atoms, composition, parent, num_zones=10):
    """Internal helper for shuffling metal lattice while protecting interstitials."""
    # 1. Define what we ARE allowed to shuffle (Metal species only)
    allowed_to_shuffle = list(composition.keys()) # Exclude interstitials
    
    y_pos = [atom.position[1] for atom in atoms]
    min_y, max_y = min(y_pos), max(y_pos)
    zone_h = (max_y - min_y) / num_zones
    zones = {i: [] for i in range(num_zones)}

    # 2. Only add an atom to a zone if it is a metal species
    # This automatically ignores 'B', 'C', etc.
    for i, atom in enumerate(atoms):
        if atom.symbol in allowed_to_shuffle:
            z_idx = min(int((atom.position[1] - min_y) // zone_h), num_zones - 1)
            zones[z_idx].append(i)

    # 3. Proceed with shuffling only the filtered indices
    for z_indices in zones.values():
        if not z_indices: continue
        
        z_symbols = []
        for s, frac in composition.items():
            z_count = int(round(len(z_indices) * frac))
            z_symbols.extend([s] * z_count)
        
        diff = len(z_indices) - len(z_symbols)
        if diff > 0:
            z_symbols.extend([parent] * diff)
        
        np.random.shuffle(z_symbols)
        
        for idx, sym in zip(z_indices, z_symbols):
            atoms[idx].symbol = sym
            
    return atoms

def initial_relax(system, calculator):
    system.calc = calculator
    system.wrap()

    logfile = 'initial_relax.log'
    if os.path.exists(logfile):
        os.remove(logfile)

    FIRE(system, logfile=logfile).run(fmax=0.01, steps=350)
    FIRE(UnitCellFilter(system), logfile=logfile).run(fmax=0.01, steps=350)
    return system, system.get_potential_energy()

def relax_config(system, move_type, calculator):
    system.calc = calculator
    system.wrap()
    # Log special moves for troubleshooting
    # remove existing logs before adding new ones
    if os.path.exists(f"{move_type}_relax.log"):
        os.remove(f"{move_type}_relax.log")

    log = f"{move_type}_relax.log" if move_type in ['insert', 'delete', 'new_host', 'flip', 'spectral'] else None

    qn = MDMin(system, logfile=log)
    qn.run(fmax=0.05, steps=30)

    dyn = FIRE(system, logfile=log)
    dyn.run(fmax=0.05, steps=50)
    return system, system.get_potential_energy()

def place_near_host(atoms, host_index, bc_index):
    """
    Simplified seeding: Places the interstitial (bc_index) in a valid 
    geometric pocket specifically around one metal host (host_index).
    """
    default_min_distance = 1.0
    relaxed_min_distance = 0.75
    
    # Pre-defined high-symmetry displacement vectors (Octahedral/Tetrahedral directions)
    displacement_distance = 1.0
    sqrt_dist = 1.2 * np.sqrt(2)

    displacement_vectors = [
        # Along principal axes
        np.array([ sqrt_dist, 0.0, 0.0 ]),
        np.array([-sqrt_dist, 0.0, 0.0 ]),
        np.array([ 0.0, sqrt_dist, 0.0 ]),
        np.array([ 0.0,-sqrt_dist, 0.0 ]),
        np.array([ 0.0, 0.0, sqrt_dist ]),
        np.array([ 0.0, 0.0,-sqrt_dist ]),
    ]

    def is_position_valid(pos, min_dist):
        """Checks if the candidate position overlaps with any atom (except host/itself)."""
        for idx, atom in enumerate(atoms):
            if idx in [host_index, bc_index]:
                continue
            # Basic distance check; for periodicity, use atoms.get_distances if needed
            if np.linalg.norm(atom.position - pos) < min_dist:
                return False
        return True

    # Get the coordinate of our target metal host
    host_pos = atoms[host_index].position
    
    # Shuffle vectors to ensure stochasticity in site selection
    random.shuffle(displacement_vectors)

    # 1. Attempt placement with standard distance (1.0 A)
    for disp in displacement_vectors:
        candidate_pos = host_pos + disp
        if is_position_valid(candidate_pos, default_min_distance):
            return candidate_pos

    # 2. Relax criteria if first pass fails (0.75 A)
    # Useful for tight Grain Boundary cores
    for disp in displacement_vectors:
        candidate_pos = host_pos + disp
        if is_position_valid(candidate_pos, relaxed_min_distance):
            return candidate_pos

    # 3. Last Resort: If even relaxed fails, return None (Move Rejected)
    return None

def insert_deletion_move_logic(system, voids_path, species, move_type):
    """Handles particle exchange logic."""
    voids = read(voids_path)
    void_coords = voids.positions.copy()
    additives = list(species.keys())

    if move_type == "insert":
        if len(void_coords) == 0: return system, 0.0, 0.0, None
        idx = random.randint(0, len(void_coords) - 1)
        symbol = random.choice(additives)
        system += Atoms(symbol, positions=[void_coords[idx]])
        new_voids = voids.copy(); del new_voids[idx]
        prefactor = len(void_coords) / sum(1 for a in system if a.symbol == symbol)
        return system, species[symbol], prefactor, new_voids

    elif move_type == "delete":
        dopant_indices = [i for i, a in enumerate(system) if a.symbol in additives]
        if not dopant_indices: return system, 0.0, 0.0, None
        idx = random.choice(dopant_indices)
        symbol, pos = system[idx].symbol, system[idx].position.copy()
        system.pop(idx)
        new_voids = voids.copy(); new_voids += Atoms("H", positions=[pos])
        prefactor = (len(dopant_indices)) / (len(void_coords) + 1)
        return system, species[symbol], prefactor, new_voids

def decorate_environment(system, solute, num_to_replace, int_idx, config):
    """
    Systematically replaces the nearest metal neighbors of an interstitial 
    with a solute species.
    """
    from core.moves import MoveSelector
    ms = MoveSelector(config)
    
    # Get neighbors within the specified spectral cutoff
    neighbors = ms._get_neighbors(system, int_idx, multiplier=1.2)
    
    # Sort neighbors by distance to the interstitial and take the nearest 6
    # (Note: _get_neighbors using natural_cutoffs generally returns nearest shells)
    metal_neighbors = [n for n in neighbors if system[n].symbol not in ms.interstitials]
    
    # Decorate the closest neighbors up to num_to_replace
    selected = metal_neighbors[:min(num_to_replace, len(metal_neighbors))]
    
    for idx in selected:
        system[idx].symbol = solute
        
    return system

def get_effective_mu(base_mu, c_target, temperature):
    """Internalizes the logic formerly in execute.py"""
    if c_target is None: return base_mu
    ln_term = np.log(c_target / (1.0 - c_target))
    return base_mu + (K_B * temperature * ln_term)

def run_md(atoms, calculator, temp, steps):
    """
    Standard NVT Langevin MD routine.
    """
    atoms.calc = calculator
    atoms.wrap()
    
    # Standard Hybrid MD parameters
    dt = 1.0 * units.fs
    friction = 0.02
    
    # erase existing logs before adding new ones
    if os.path.exists('nvt_md.log'):
        os.remove('nvt_md.log')
        
    dyn = Langevin(atoms, timestep=dt, temperature_K=temp, 
                   friction=friction, logfile='nvt_md.log', 
                   loginterval=50)
    dyn.run(steps)
    
    return atoms, atoms.get_potential_energy()


from ase.mep import NEB, NEBTools
import matplotlib.pyplot as plt
import logging

def run_neb_calculation(atoms_0, atoms_1, calculator_setup_func, symbol):
    """
    NEB wrapper using a 3-point guided interpolation (Start -> Saddle -> End).
    Dynamically identifies the interstitial atom by its symbol.
    """

    # Define the log path
    neb_log = 'neb_current.log'

    # 1. Delete the log file if it exists to start fresh for this event
    if os.path.exists(neb_log):
        os.remove(neb_log)

    # 1. Dynamically find the index of the interstitial atom
    # This is more robust than assuming index 0
    indices = [atom.index for atom in atoms_0 if atom.symbol == symbol]
    
    if len(indices) == 0:
        logging.error(f"NEB Error: Interstitial symbol '{symbol}' not found in system.")
        return 0.0, 0.0
    if len(indices) > 1:
        logging.warning(f"NEB Warning: Multiple '{symbol}' atoms found. Using first instance (index {indices[0]}).")
    
    idx = indices[0]

    # 2. Setup 8 images (initial, 6 intermediate, final)
    images = [atoms_0.copy()]
    for _ in range(6):
        images.append(atoms_0.copy())
    images.append(atoms_1.copy())

    # 3. Linear Interpolation of the Interstitial Only
    pos0 = atoms_0.positions[idx]
    pos1 = atoms_1.positions[idx]
    
    n_images = len(images)
    for i in range(1, n_images - 1):
        # Standard linear interpolation: r(i) = r0 + fraction * (rf - r0)
        fraction = i / (n_images - 1)
        interstitial_pos = (1 - fraction) * pos0 + fraction * pos1
        images[i].positions[idx] = interstitial_pos

    # 4. Initialize NEB and Calculators
    neb = NEB(images, k=0.5, climb=True, parallel=True)
    
    # Assign unique calculators to each image to prevent state conflicts
    for img in images:
        img.calc = calculator_setup_func()

    # 5. Run Optimization
    qn = MDMin(neb, dt=0.02, logfile=neb_log)
    qn.run(fmax=0.05, steps=50)

    qn = FIRE(neb, maxstep=0.01, logfile=neb_log)
    qn.run(fmax=0.05, steps=300)

    # 6. Extract forward/reverse barriers
    nebtools = NEBTools(images)
    dEf, _ = nebtools.get_barrier()

    E0 = images[0].get_potential_energy()
    Ef = images[-1].get_potential_energy()

    dEr = dEf - (Ef - E0)

    fig = nebtools.plot_band()
    fig.savefig(f'data/neb/MEP.png', bbox_inches='tight')
    plt.close(fig)

    return dEf, dEr