# Example: Standard Hybrid MC/MD Simulation

This example demonstrates a standard "equilibration" run for a Ni-Cr alloy system. It utilizes the **Hybrid Monte Carlo/Molecular Dynamics** approach to sample chemical ordering while simultaneously allowing the lattice to relax mechanically.



## ⚙️ Simulation Setup

The provided `config.yaml` is configured for a **Canonical Ensemble** run with the following characteristics:

### 1. Hybrid Routine Control
* **Monte Carlo (MC)**: The simulation performs 35,000 MC steps to sample atomic swaps.
* **Molecular Dynamics (MD)**: Every 2,000 MC steps, the routine pauses to run 20,000 steps of MD relaxation. This ensures that the lattice geometry stays consistent with the changing chemical environment.
* **Calculator**: Powered by the **PFP** machine-learning potential.

### 2. System Construction
* **Cell Type**: Automatically generates an **FCC** supercell ($8 \times 8 \times 8$ unit cells).
* **Composition**: Initialized with **85% Ni** and **15% Cr**.
* **Randomization**: Atomic symbols are scrambled at the start to ensure a non-biased initial state.

### 3. Active Ensembles
* **Canonical**: Enabled to allow for position-based swaps.
* **Semi-Grand & Grand Canonical**: Disabled for this specific example.

---

## 💡 Usage Notes

### Switching to Pure Monte Carlo
If you do not require MD relaxations and wish to perform a **Pure MC** run (e.g., for faster sampling of chemical ordering on a rigid lattice), modify the `simulation` block as follows:

```yaml
simulation:
  hybrid_md:
    enabled: False  # Disables the MD thermostat/relaxation calls