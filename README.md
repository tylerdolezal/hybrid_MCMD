# Hybrid MC/MD Simulation & Spectral Mapping Tool

This framework is designed for exploring the thermodynamics and kinetics of interstitials in complex environments (e.g., grain boundaries, surfaces, and bulk alloys). It integrates stochastic Monte Carlo (MC) sampling with molecular dynamics (MD) relaxations and systematic energy mapping using high-fidelity potentials like **PFP**, **CHGNet**, or **EAM**.

## 🚀 Core Functionalities

### 1. Hybrid MC/MD Simulation
Performs random-walk sampling of atomic configurations to reach thermodynamic equilibrium:
* **Canonical Ensemble**: Position swaps for metal atoms and host-to-host diffusion for interstitials.
* **Semi-Grand Canonical**: Metal species identity flips (e.g., $Ni \leftrightarrow Cr$) based on chemical potential differences ($\Delta\mu$).
* **Grand Canonical**: Stochastic insertion and deletion of interstitial species (e.g., B, C, N, O) based on target concentrations ($c_{target}$) and reservoir potentials.

### 2. Spectral Mapping & NEB
A systematic mode for mapping the energy landscape of a single interstitial through a specific diffusion network:
* **Random Walk Diffusion**: Moves an interstitial between neighboring sites identified in a `voids_file`.
* **Deterministic Decoration**: Systematically replaces the nearest metal neighbors of the interstitial with a solute (e.g., $Cr$) to measure alloying effects.
* **Kinetic Barriers**: Automatically calculates forward ($dE_f$) and reverse ($dE_r$) activation barriers using the Nudged Elastic Band (NEB) method.

---

## 🛠 Directory Structure

* `main.py`: The primary entry point. Handles mode dispatching.
* `core/simulation.py`: Contains the logic for `HybridSimulation` and the `SpectralCollector` subclass.
* `core/thermo.py`: Logic for relaxations, NEB image generation, and local environment decoration.
* `core/moves.py`: Logic for executing the various MC moves. 
* `drivers/`: Potential-specific wrappers (e.g., `pfp_driver.py`) that return ASE-compatible calculators.
* `utils/`: Configuration parsing, standard I/O, and dual-output logging.
* `data/`: Output directory for MC-MD data, `spectral_log.csv`, and NEB trajectory data.

---

## 📝 Configuration (`config.yaml`)

The tool uses a nested YAML structure. For **Spectral Mode**, ensure the following parameters are set:

```yaml
simulation:
  num_mc_steps: 35000
  temperature: 1200        # Global temperature for MC acceptance and MD thermostat
  snapshot_every: 200
  continue_run: False
  potential_style: "pfp" # Options: pfp, chgnet, eam
  
  # Explicit control over the MD portion of the hybrid routine
  hybrid_md:
    enabled: False         
    interval: 500          # MC steps between MD calls
    steps: 200              # Number of MD steps per call
    
spectral:
  enabled: True            # Toggle systematic mapping vs. random simulation
  do_neb: True             # Enable activation barrier calculations
  jump_cutoff: 2.75         # Max distance (Å) for a valid diffusion hop
  solute: "Cr"             # Solute used for deterministic decoration shell

system:
  use_custom_cell: True    # Required: Loads 'POSCAR-custom'

ensembles:
  grand:
    enabled: True          # Must be True to provide species and void data
    voids_file: "gb_voids-POSCAR-S5" 
    additives:
      B: { base_mu: -6.46 } # Species to be mapped