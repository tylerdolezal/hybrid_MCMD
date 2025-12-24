# Void Mapping & Connectivity Analysis: $\Sigma 5$ Symmetric Tilt Grain Boundary

This example demonstrates how to pre-process a complex atomic structure—specifically a **pure Nickel $\Sigma 5$ (210) symmetric tilt grain boundary (STGB)**—to identify valid interstitial diffusion sites. The algorithm uses Voronoi tessellation and hierarchical clustering to map the "diffusion network" required for Spectral Mapping and NEB calculations.

![alt text](image.png)

## 🔬 Algorithm Overview

The mapping process follows four distinct stages to transform raw atomic coordinates into a clean network of interstitial nodes:

### 1. Spatial Bounding Box Filtering
The script reads a `POSCAR-S5-gb` file where the grain boundary interface is centered along the **y-axis**. 
* **Targeting**: Using the `S5_bounds` dictionary, the algorithm restricts its search to a specific "slice" of the simulation cell.
* **Focus**: This allows the user to target the disordered GB core specifically (e.g., $y \in [29.0, 34.0]$) while ignoring bulk regions where interstitial sites are already well-known.



### 2. Voronoi Vertex Detection & Radius Filtering
The algorithm detects low-density regions using **Voronoi Tessellation**:
* **Vertices as Candidates**: Every Voronoi vertex is a point equidistant from at least four atoms, making them the geometric centers of potential octahedral or tetrahedral holes.
* **Radius Threshold (`rmin`)**: Only vertices that are at least $1.54$ Å (for this Ni example) from the nearest atom are kept. This effectively filters out small, unstable "tight" spaces and focuses on viable interstitial pockets.

### 3. Hierarchical COM Collapse
In distorted GB regions, multiple Voronoi vertices often cluster around a single physical void.
* **Agglomerative Clustering**: The algorithm uses a fixed distance threshold (default $0.8$ Å) to group these sub-vertices.
* **Center of Mass (COM)**: Each cluster is collapsed into a single representative "Void Marker" located at the cluster's center of mass. This prevents "double-counting" sites during Monte Carlo moves.

### 4. Network Connectivity Analysis
The final step ensures the detected voids form a continuous path for diffusion.
* **Natural Shell Detection**: The script finds the average nearest-neighbor distance between voids.
* **Dynamic Buffer**: It adds a user-defined buffer (default $0.5$ Å) to this average to find the "natural" connectivity shell of the network.



---

## 📤 Output Files

Upon execution, the script generates two primary VASP-formatted files:

1.  **`voids-POSCAR-S5`**: Contains *only* the identified void markers (represented as `H` atoms). This file is used as the `voids_file` input in your `config.yaml`.
2.  **`full_voids-POSCAR-S5`**: A visualization file containing both the original metal lattice and the new void markers. Open this in **OVITO** or **VESTA** to verify that the voids are correctly positioned within the grain boundary interface.

---

## 📏 Determining the `jump_cutoff`

The script provides a console report to help you accurately configure the `spectral` block in `config.yaml`. 

**The Connectivity Report looks like this:**
```text
=============================================
      VOID NETWORK CONNECTIVITY REPORT
=============================================
Natural NN Average:      2.485 Å
Search Threshold:        2.985 Å
Average Edge Distance:   2.512 Å
MAX EDGE (Suggested):    2.715 Å
---------------------------------------------
RECOMMENDED jump_cutoff: 2.80 Å
=============================================