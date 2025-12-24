import numpy as np
from ase import Atoms
from ase.atom import Atom
import matplotlib.pyplot as plt

def insert_fcc_interstitials(
    atoms,
    interstitial_element='H',
    min_dist=1.0,
    bounds=None
):
    """
    Places interstitial atoms at ideal FCC octahedral and tetrahedral sites
    throughout the given atomic cell, restricted to a bounding box if provided.

    Parameters:
        atoms: ASE Atoms object
        interstitial_element: Element symbol for interstitials (e.g., 'H')
        min_dist: Minimum allowed distance to existing atoms (Å)
        bounds: Optional spatial bounding box:
            {
                'xlo': float, 'xhi': float,
                'ylo': float, 'yhi': float,
                'zlo': float, 'zhi': float
            }

    Returns:
        ASE Atoms object with interstitials added at valid sites
    """

    cell = atoms.get_cell()
    lengths = np.linalg.norm(cell, axis=1)
    lattice_constant = 3.524  # approximate FCC spacing
    nx, ny, nz = [int(np.ceil(length / lattice_constant)) for length in lengths]

    # FCC ideal interstitial site types in fractional coords
    octa_sites = [
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
        [0.5, 0.5, 0.5],
    ]
    tetra_sites = [
        [0.25, 0.25, 0.25],
        [0.75, 0.25, 0.25],
        [0.25, 0.75, 0.25],
        [0.25, 0.25, 0.75],
    ]
    
    interstitials = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for frac, symbol in zip(octa_sites + tetra_sites, ['C']*len(octa_sites) + ['B']*len(tetra_sites)):
                    frac_pos = np.array([i + frac[0], j + frac[1], k + frac[2]])
                    cart_pos = np.dot(frac_pos, cell)/10 # dot product of fractional vector with cell matrix

                    # Bounding box check
                    if bounds:
                        x, y, z = cart_pos
                        if not (bounds['xlo'] <= x <= bounds['xhi'] and
                                bounds['ylo'] <= y <= bounds['yhi'] and
                                bounds['zlo'] <= z <= bounds['zhi']):
                            continue

                    # Distance to existing atoms
                    if all(np.linalg.norm(atom.position - cart_pos) > min_dist for atom in atoms):
                        interstitials.append(Atom(symbol, cart_pos))

    # Combine
    new_atoms = atoms.copy()
    for interstitial in interstitials:
        new_atoms.append(interstitial)

    return new_atoms


from scipy.spatial import Voronoi, cKDTree
from sklearn.cluster import AgglomerativeClustering

def collapse_voids_to_COM_hierarchical(void_atoms, cutoff=0.8, element='H'):
    """
    Collapses tightly packed void atoms into center-of-mass (COM) representatives
    using hierarchical clustering with a fixed cutoff.

    Parameters:
        void_atoms: ASE Atoms object with void markers (e.g., 'X')
        cutoff: max distance (Å) for grouping atoms into a single void cluster
        element: element symbol for resulting COM atoms

    Returns:
        ASE Atoms object with one atom per void site (COM marker)
    """
    positions = void_atoms.get_positions()
    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=cutoff,
        linkage='single'
    ).fit(positions)

    labels = clustering.labels_
    collapsed_atoms = []

    for label in np.unique(labels):
        cluster_pos = positions[labels == label]
        com = cluster_pos.mean(axis=0)
        collapsed_atoms.append(Atom(element, com))

    return Atoms(collapsed_atoms)


def find_voids_from_atoms(atoms, bounds, min_radius=1.0, element='H'):
    """
    Detect voids (low-density regions) using Voronoi vertices within a bounding box.

    Parameters:
        atoms: ASE Atoms object
        bounds: dict with Cartesian bounds {'xlo','xhi','ylo','yhi','zlo','zhi'}
        min_radius: Minimum allowed radius to nearest atom (void size threshold)
        return_as_atoms: if True, return as ASE Atoms object; else return positions
        element: element symbol to assign to void markers (e.g., 'X')

    Returns:
        Atoms object with 'X' atoms at void centers, or list of positions
    """
    positions = atoms.get_positions()
    vor = Voronoi(positions)
    tree = cKDTree(positions)

    void_atoms = []
    for vertex in vor.vertices:
        x, y, z = vertex
        # Bounding box filter
        if not (bounds['xlo'] <= x <= bounds['xhi'] and
                bounds['ylo'] <= y <= bounds['yhi'] and
                bounds['zlo'] <= z <= bounds['zhi']):
            continue

        # Distance to nearest atom
        # Filter by distance and also "connectivity" to other voids
        dist, _ = tree.query(vertex)
        if dist > min_radius:
            void_atoms.append(Atom(element, vertex))


    print(len(void_atoms))
    # Return combined structure
    new_atoms = atoms.copy()

    for v in void_atoms:
        new_atoms.append(v)
    
    # Extract just the voids ("H" atoms)
    voids_only = new_atoms[[atom.symbol in ['B', 'C', 'H'] for atom in new_atoms]]
    
    # Collapse void clusters into COM markers
    collapsed_voids = collapse_voids_to_COM_hierarchical(voids_only)
    print(len(collapsed_voids))
    # Return combined structure
    new_atoms = atoms.copy()

    for v in collapsed_voids:
        new_atoms.append(v)
    
    # Extract just the voids ("H" atoms)
    voids_only = new_atoms[[atom.symbol in ['B', 'C', 'H'] for atom in new_atoms]]

    # Get void positions
    void_positions = np.array([atom.position for atom in voids_only])

    # Build KDTree and query neighbors within 2.0 Å
    tree = cKDTree(void_positions)
    pairs = tree.query_ball_tree(tree, r=3.0)

    # Keep only voids that have at least one *other* void nearby
    valid_indices = [i for i, neighbors in enumerate(pairs) if len(neighbors) > 1]

    # Filter voids_only down to those with nearby neighbors
    voids_only = voids_only[valid_indices]

    print(len(voids_only))

    return voids_only, new_atoms

def analyze_void_connectivity(voids_only, buffer=0.5):
    """
    Analyzes distances between neighboring voids to determine the 
    optimal 'jump_cutoff' for spectral mapping.
    
    Parameters:
        voids_only: ASE Atoms object containing void markers.
        buffer: Extra Ångströms to add to the average nearest-neighbor distance.
    """
    positions = voids_only.get_positions()
    if len(positions) < 2:
        print("Insufficient voids to analyze connectivity.")
        return None
        
    tree = cKDTree(positions)
    
    # 1. Determine the "natural" first-neighbor shell distance
    # Query k=2 because k=1 is the atom itself (distance 0)
    dists, _ = tree.query(positions, k=2)
    nn1_distances = dists[:, 1]
    
    natural_avg = np.mean(nn1_distances)
    dynamic_threshold = natural_avg + buffer
    
    # 2. Find all edges within this threshold
    pairs = tree.query_ball_tree(tree, r=dynamic_threshold)
    
    edge_distances = []
    for i, neighbors in enumerate(pairs):
        for neighbor_idx in neighbors:
            if i < neighbor_idx:  # Prevent double-counting and self-pairing
                d = np.linalg.norm(positions[i] - positions[neighbor_idx])
                edge_distances.append(d)
    
    if not edge_distances:
        print(f"No connections found within a {dynamic_threshold:.2f} Å threshold.")
        return None

    avg_edge = np.mean(edge_distances)
    max_edge = np.max(edge_distances)
    
    print("\n" + "="*45)
    print("      VOID NETWORK CONNECTIVITY REPORT")
    print("="*45)
    print(f"Total Voids Found:       {len(positions)}")
    print(f"Minimum NN Distance:     {np.min(nn1_distances):.3f} Å")
    print(f"Natural NN Average:      {natural_avg:.3f} Å")
    print(f"Search Threshold:        {dynamic_threshold:.3f} Å")
    print(f"Average Edge Distance:   {avg_edge:.3f} Å")
    print(f"MAX EDGE (Suggested):    {max_edge:.3f} Å")
    print("-" * 45)
    print(f"RECOMMENDED jump_cutoff: {np.ceil(max_edge * 10)/10:.2f} Å")
    print("="*45 + "\n")
    
    return max_edge

from ase.io import read, write

bulk_bounds = {
                'xlo': 5.0, 'xhi': 19.0,
                'ylo': 5.0, 'yhi': 19.0,
                'zlo': 5.0, 'zhi': 19.0
            } # 'ylo': 26.5, 'yhi': 36 ; generates a few bulk voids


atoms = read("POSCAR-bulk (pure Ni)", format="vasp")

rmin = 1.54 # this gets rid of tetrahedral sites in the bulk
voids_only, full_atoms = find_voids_from_atoms(atoms, bulk_bounds, min_radius=rmin)
write("voids-POSCAR-bulk", voids_only, format="vasp", direct=True, sort=True)
write("full_voids-POSCAR-bulk", full_atoms, format="vasp", direct=True, sort=True)

analyze_void_connectivity(voids_only)