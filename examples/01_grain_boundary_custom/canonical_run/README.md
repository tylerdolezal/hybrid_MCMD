# Example: Canonical MC/MD with Custom Cell & Interstitials

This example demonstrates a research-grade simulation using a **custom geometric input** (e.g., a Grain Boundary or Surface) while maintaining a specific chemical composition through Monte Carlo swaps.


## ⚙️ Simulation Setup

The configuration is optimized for a **Canonical Ensemble** run using a pre-defined lattice structure:

### 1. Custom Geometry Logic
* **`use_custom_cell: True`**: The simulation bypasses internal crystal generation and directly loads **`POSCAR-custom`**.
* **Interstitial Protection**: While the metal lattice is initialized to **85% Ni and 15% Cr**, the randomization scheme is configured to ignore existing interstitial species (e.g., Boron). This ensures your manually placed interstitial remains in its starting site while the surrounding alloy environment equilibrates.

### 2. Hybrid Routine Control
* **Monte Carlo (MC)**: Performs 35,000 steps to sample the chemical ordering of Ni and Cr.
* **Molecular Dynamics (MD)**: Runs 20,000 steps of MD relaxation every 2,000 MC steps. This is critical for custom cells to allow the grain boundary or defect structure to relax mechanically as the local chemistry changes.
* **Temperature**: Constant at **600 K** for both the Metropolis criterion and the MD thermostat.

### 3. Active Ensembles
* **Canonical**: Enabled for position-based metal-metal swaps and interstitial host-to-host jumps.
* **Semi-Grand & Grand Canonical**: Disabled to maintain a fixed number of atoms and species.

---

## 💡 Usage Notes

### Maintaining Composition
The `composition` block is used to decorate the `POSCAR-custom` file upon loading. If the initial file is pure Ni, the code will automatically flip the required number of sites to Cr to reach the 15% target before starting the MC loop.

### Pure Monte Carlo (Static Lattice)
If you wish to keep the `POSCAR-custom` coordinates strictly fixed and only sample chemical ordering, disable the MD portion:

```yaml
simulation:
  hybrid_md:
    enabled: False  # Disables lattice relaxation