# Example: Spectral Mapping & NEB (Grain Boundary Complement)

This example demonstrates the **Spectral Mapping** mode applied specifically to a **$\Sigma 5$ Grain Boundary (GB)** environment. It is designed to map the energetic and kinetic landscape of an interstitial moving along the disordered interface of a grain boundary.



## 🔬 Spectral Logic: Grain Boundary Application

In the grain boundary complement, the systematic survey focuses on how the unique geometry of the interface—combined with localized solute decoration—alters the diffusion path of an interstitial.

### 🛠️ Configuration Adjustments for GB
* **`use_custom_cell: True`**: This mode is essential for loading the specific $\Sigma 5$ GB geometry from `POSCAR-custom`.
* **`voids_file: "voids-POSCAR-S5"`**: The simulation utilizes the specific void network generated for the $\Sigma 5$ interface to ensure jumps are restricted to the GB core.
* **`jump_cutoff: 2.0`**: This value has been optimized to **2.0 Å** based on the recommendation from the `void_generation` script. This ensures the "Spectral Sampler" identifies the high-density hopping paths unique to the GB structure without skipping valid neighboring sites.
* **Ensemble Lock**: As with the bulk version, all ensembles except **`canonical`** remain disabled to ensure a controlled, non-stochastic survey of the energy landscape.

### 🚶 The Mapping Process
1. **Iterative Decoration**: The sampler systematically "decorates" the nearest-neighbor shell of the interstitial with the chosen solute (e.g., `Cr`) from **0 to 6 neighbors** at every jump site.
2. **Kinetic Barriers (`do_neb: True`)**: The tool calculates Climbing Image Nudged Elastic Band (CI-NEB) barriers for every jump in the GB network. This reveals how grain boundary strain and local chemistry synergistically affect the activation energy for diffusion.
3. **Static Spectra (`do_neb: False`)**: If NEB is toggled off, the tool simply records the total energy at each site, providing a "segregation spectrum" that shows where the interstitial is most stable along the interface.



---

## 📏 Parameters & Execution

* **`num_mc_steps: 100`**: Defines the number of successful jumps within the GB void network.
* **`snapshot_every: 1`**: Records every unique configuration to capture the full spectrum of local environments encountered.
* **Lattice Dynamics**: **MD is disabled** (`hybrid_md: enabled: False`) to maintain the integrity of the spectral mapping on a static, controlled lattice.

---

## 💡 Usage Notes

### Analyzing GB Trapping
By comparing this data to the bulk FCC example, you can quantify the "trapping" strength of the grain boundary. If the energy barriers ($dE_f$) are significantly higher or the total energies significantly lower in the GB core, the interface acts as a kinetic or thermodynamic sink for the interstitial.

### Running the Example
```bash
python main.py --config config.yaml