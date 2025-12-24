# Example: Fully Active Co-Doped GB Simulation (GC + SGC + Canonical)

This example represents the most advanced "all-active" workflow in the framework. It simulates a grain boundary (GB) in a multi-component environment where metal segregation ($\Delta\mu_{Ni-Cr}$) and interstitial co-doping (Boron and Carbon) occur simultaneously under full thermodynamic equilibrium.


## 🔬 Simulation Logic

This setup models a "Grand Reservoir" where the grain boundary can exchange metal identities with a bulk alloy reservoir and exchange interstitials with an external chemical atmosphere.

### 1. Interstitial Thermodynamics (GC)
* **Species**: Trialing the insertion and deletion of both **Boron** and **Carbon**.
* **Reference State**: The `base_mu` values (-6.46 eV for B, -7.51 eV for C) are derived from single-atom insertion energies into a **pure Ni FCC bulk octahedral site**. This provides a stable baseline that remains consistent even as the global Chromium concentration varies.
* **Fermi-Dirac Modulation**: The code uses a target concentration (`c_target`) and simulation temperature to calculate an **Effective Chemical Potential** ($\mu_{eff}$), accounting for site occupancy statistics:
  $$\mu_{eff} = \mu_{base} + k_B T \ln\left(\frac{c_{target}}{1 - c_{target}}\right)$$

### 2. Metal Segregation (SGC)
* **Converged $\Delta\mu$**: Utilizes a pre-calibrated `delta_mu: 0.936` to maintain a global concentration of **25 at% Cr** at **300 K**. 
* **Dynamic Exchange**: While the global reservoir is fixed, Cr is free to preferentially segregate to the grain boundary core based on the machine-learning potential (PFP) energy landscape.

### 3. Structural Relaxation (Hybrid MD)
* **Lattice Breathing**: `hybrid_md` is fully engaged, running 20,000 steps of MD relaxation every 2,000 MC steps. This is essential for co-doped systems to allow the GB to expand or contract locally as interstitials are inserted into void sites.

---

## ⚙️ Simulation Setup

* **Custom Cell**: Loads `POSCAR-custom` (the $\Sigma 5$ GB interface).
* **Void Network**: Uses `voids-POSCAR-S5` to define the legal trial sites for B and C insertions.
* **Ensembles**:
    * **Canonical**: Handles metal swaps and long-range interstitial diffusion.
    * **Semi-Grand**: Handles $Ni \leftrightarrow Cr$ identity transmutations.
    * **Grand**: Handles the stochastic loading of B and C.


---

## 💡 Usage Notes

### ⚠️ Potential & Composition Dependency
Be advised that while the `base_mu` is derived from pure Ni to remain "transferable" across different alloy concentrations, the resulting occupancy will still vary if you change the interatomic potential (e.g., switching from PFP to CHGNet) or the temperature, as these affect the $k_B T$ modulation and the underlying PES.

### Analyzing Co-Doping Synergies
With both B and C active, you can observe **competitive occupancy** or **synergistic stabilization**. Does the presence of segregated Cr at the grain boundary make the insertion of Boron more or less energetically favorable?

### Running the Example
```bash
python main.py --config config.yaml