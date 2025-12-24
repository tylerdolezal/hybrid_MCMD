# Configuration Guide: Interstitial Mapping & Hybrid MC/MD

This guide provides a detailed breakdown of the parameters within `config.yaml`. This file controls the thermodynamic ensembles, kinetic mapping, and potential energy surface (PES) evaluation settings.

---

## 🟢 Simulation Settings
These parameters define the global execution environment.

* **`num_mc_steps`**: The total number of Monte Carlo trials. In **Spectral Mode**, this acts as the maximum number of successful diffusion "jumps" to record.
* **`temperature`**: Set to **1200 K** by default to facilitate crossing of higher energy barriers in standard MC. For **Spectral Mode**, this temperature is used by the thermo module to calculate the effective chemical potential ($\mu_{eff}$) based on target concentrations.
* **`potential_style`**: Determines the backend driver (e.g., **pfp**, **chgnet**, **eam**). This must match the installed API clients in your environment.
* **`hybrid_md`**: 
    * **`enabled`**: Should be `False` for Spectral Mapping to ensure static energy evaluations. 
    * **`interval/steps`**: Define how often the system undergoes molecular dynamics to relax the global lattice during a standard MC run.

---

## 🔵 System Construction
Defines how the initial atomic geometry is handled.

* **`use_custom_cell`**: When `True`, the code ignores the `crystal_shape` and `supercell_size` parameters and directly loads **`POSCAR-custom`**. This is required for Grain Boundary (GB) studies.
* **`composition`**: The target ratio of metal atoms (e.g., `{Ni: 1.0, Cr: 0.0}`). In Spectral Mode, this is ignored in favor of the systematic decoration logic.
* **`randomize_initial`**: If `True`, the metal symbols are shuffled at step 0. Keep this `False` for Spectral scans where you want a specific starting configuration.

---

## 🟡 Ensembles
The framework supports three thermodynamic ensembles that can be toggled independently.

### 1. Canonical
Handles the exchange of atom positions (swaps) and diffusion. In Spectral Mode, this logic is bypassed by the `SpectralCollector`.

### 2. Semi-Grand Canonical (SGC)
* **`mode`**: Options like `fixed_mu` determine how species identity is swapped.
* **`metal_library`**: The elements available for identity flips (e.g., Ni and Cr).
* **`delta_mu`**: The chemical potential difference between species, provided as a file or float.

### 3. Grand Canonical (GC)
Used for managing interstitial species.
* **`voids_file`**: Path to the file containing valid interstitial sites (Octahedral/Tetrahedral pockets).
* **`additives`**: 
    * **`base_mu`**: The reference chemical potential for the species (e.g., $-6.46$ eV for Boron). 
    * **`c_target`**: The target concentration used to calculate effective $\mu$ at a given temperature.

---

## 🟣 Spectral Mapping (V2 Workflow)
This block controls the systematic exploration of the grain boundary energy landscape.

* **`enabled`**: If `True`, the simulation ignores random MC moves and begins a systematic walk through the `voids_file`.
* **`do_neb`**: Activates the Nudged Elastic Band calculation for every successful jump between sites.
* **`jump_cutoff`**: Set to **2.75 Å**. This is the maximum distance allowed for a "neighboring" void to be considered a valid jump target.
* **`solute`**: The species used to "decorate" the local environment (e.g., **Cr**).



---

## 📂 Example: Mapping Boron in Ni-Cr
To perform a spectral scan of Boron at a Grain Boundary:
1. Ensure `spectral: enabled` is `True`.
2. Set `use_custom_cell` to `True` and provide your `POSCAR-custom`.
3. Provide a `voids_file` corresponding to that specific POSCAR.
4. Set your solute to `Cr` and ensure `grand: enabled` is `True` to load the `B` parameters.