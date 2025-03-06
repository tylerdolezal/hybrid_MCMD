# Hybrid MCMD Input File Tutorial

This tutorial explains the structure and usage of the `input_file` for the **Hybrid Monte Carlo Molecular Dynamics (MCMD) Routine**.

## Overview
The input file contains key parameters that define the simulation settings, including composition, crystal structure, Monte Carlo (MC) and Molecular Dynamics (MD) parameters, and additional modifiers like grain boundaries, vacancies, and additives.

Each line in the file follows a `key: value` format, where:
- Keys define the parameter type.
- Values specify the corresponding settings.

## Parameter Descriptions

### Composition
```plaintext
composition: Ni=0.7, Cr=0.3
```
Defines the atomic fraction of elements in the system. In this case, the system contains **70% Ni and 30% Cr**.

### Crystal Shape
```plaintext
crystal_shape: fcc
```
Defines the crystal structure of the simulation cell. Options include **fcc (Face-Centered Cubic), bcc (Body-Centered Cubic), and hcp (Hexagonal Close-Packed)**.

### Grain Boundary
```plaintext
grain_boundary: False
```
Specifies whether a grain boundary is included in the simulation. Set to `True` to enable grain boundary modeling.

### Randomization
```plaintext
randomize: False
```
Determines whether the initial atomic positions are randomized.

### MD Parameters
```plaintext
md_params: 0, 1073, chgnet
```
Defines the **MD execution settings**:
- `0` – Initial MD step.
- `1073` – Target temperature in Kelvin.
- `chgnet` – Interatomic potential used for MD calculations (e.g., `chgnet`, `eam`, `meam`).

### Monte Carlo Parameters
```plaintext
num_mc_steps: 6000
```
Sets the number of Monte Carlo steps to be performed.

```plaintext
md_interval: 6001
```
Determines how often MD simulations are executed within the hybrid routine.

### Simulation Cell Size
```plaintext
size: 6
```
Defines the size of the supercell in terms of unit cell repetitions.

### Execution Command
```plaintext
supcomp_command: mpirun -np 8 lmp_mpi
```
Defines the command to execute the LAMMPS simulation, specifying the number of processors (`-np 8`).

### Continue Previous Run
```plaintext
continue_run: False
```
Indicates whether the simulation should resume from a previous run.

### Additives (Dopants)
```plaintext
additives: None #[(Cr, B, 9)]
```
Defines additional elements or dopants to be introduced into the system. Here, `Cr-B` dopant with **9 atomic percent** is commented out.

### Vacancies
```plaintext
vacancies: False
```
Specifies whether atomic vacancies should be introduced in the simulation.

### Metal Library
```plaintext
metal_library: None #["Cr", "Fe", "Mo", "Nb", "Ni", "Ti"]
```
Defines the set of metallic elements available for swapping or selection.

### Surface Properties
```plaintext
surface: None #["O", "hexagonal", 7.0]
```
Defines surface conditions. Here, the `O` (oxygen) surface with a **hexagonal shape** and a thickness of **7.0 units** is commented out.

## Usage Notes
- Ensure proper syntax, as incorrect formatting may lead to parsing errors.
- Commented lines can be uncommented (`#` removed) to enable specific features.
- Default values should be reviewed before running simulations to match the intended study parameters.

---
This guide provides a structured explanation of the `input_file` format and its parameters. For any modifications, ensure consistency with the expected input syntax.
