# Hybrid Monte Carlo Molecular Dynamics (hMCMD) Routine

## Overview
The **Hybrid MCMD Routine** is designed to run Monte Carlo (MC) and Molecular Dynamics (MD) simulations in a hybrid fashion to explore material behaviors under various doping conditions. The routine allows for automated execution, statistical data collection, and efficient handling of simulation parameters. As a reminder the catalogue of MC moves are:

- **Swap atomic positions of nearest neighbor metallic constituents:** Select a pair of metallic atoms from within the nearest neighbor shell and attempt to swap their positions in the lattice.
- **Relocate a light interstitial atom near a new metallic host:** Select a light interstitial atom (B, C, H, or N) and attempt to place it near a different metallic atom within its nearest neighbor shell.
- **Introduce proximity between two light interstitial atoms:** Select two light interstitial atoms and attempt to place one of them within the first nearest neighbor shell of the other. If the first shell is fully occupied, the placement is attempted in the second nearest neighbor shell.
- **Separate two neighboring light interstitial atoms:** If such a pair exists, select a light interstitial atom and attempt to move it away from its nearest light interstitial neighbor. This move balances the proximity move, ensuring there is a chance to promote both segregation and aggregation among the light interstitials.
- **Swap a metal neighbor of a light interstitial atom:** For a selected light interstitial atom, identify one of its nearest metallic neighbors and swap that metal’s position with one of the metallic atom's own nearest neighbors.

## File Descriptions

### Main Execution Files
- **`input_file`** – User-defined settings for the simulation, including temperature, atomic species, dopant concentrations, and MC/MD parameters.
- **`execute.py`** – A script to iterate through multiple dopant concentrations or run multiple undoped simulations for statistical analysis.
- **`(p) execute.py`** – A parallelized version of `execute.py` for efficient multi-run execution.
- **`hybrid_mcmd.py`** – The core algorithm implementing hybrid Monte Carlo and Molecular Dynamics execution.

### Supporting Source Files
- **`potential.mod`** – Essential file automatically updated by the routine based on the alphabetized species found in the POSCAR files.
- **`structure.mod`** – Reads in POSCAR files for LAMMPS simulations.
- **`hybrid_md_routine.py`** – Runs MD simulations if the hybrid setting is enabled.
- **`lammps_functions.py`** – Automates the creation and modification of LAMMPS input files.
- **`myfuncs.py`** – A collection of functions used within `hybrid_mcmd.py` to streamline computations.
- **`no_MD/`** – Directory for pure Monte Carlo runs (no MD involved).
- **`w_MD/`** – Directory for hybrid Monte Carlo-Molecular Dynamics runs.
- **`chgnet_src/`** – Contains an edited version of `myfuncs.py` modified to run with CHGNet and ASE instead of LAMMPS.

## Usage Instructions
1. Modify `input_file` to specify the simulation parameters.
2. Run `execute.py` or `(p) execute.py` depending on whether a sequential or parallel execution is needed to iterate over a range of compositions or multiple simulation runs.
3. Results will be stored in the data and structures directories and in `Alloy_X{at%}_run{iteration}` if `batch_mode` is used.
4. Analyze output data for insights on material behavior.

## Dependencies
Ensure the following are installed:
- **LAMMPS**
- **ASE** 
- **Python 3**
---
This repository is actively maintained and updated to improve efficiency and accuracy in hybrid MCMD simulations. Contributions and suggestions are welcome! To help improve the code, please consider running the routine with:  

`python3 execute.py 2> error.log`

This will capture errors separately, making it easier to track and resolve any issues.
