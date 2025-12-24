# Example: Spectral Mapping & NEB (V2 Workflow)

This example demonstrates the **Spectral Mapping** mode, a systematic approach to surveying the energy landscape and kinetic barriers of an interstitial as it moves through a specific environment.


## 🔬 Spectral Logic

Unlike the stochastic Monte Carlo runs, Spectral Mapping is a deterministic survey of local environments. It systematically "decorates" the nearest-neighbor (NN) shell of an interstitial with a solute species to measure the impact of local chemistry on thermodynamics and kinetics.

### 🛠️ Configuration Requirements
* **`spectral: enabled: True`**: This redirects the `main.py` entry point to the `SpectralCollector` engine.
* **`use_custom_cell: True`**: Required. This mode must operate on a specific geometry (e.g., your 7x7x7 bulk Ni cell).
* **Ensemble Settings**: Only **`canonical`** should be enabled. While the **`voids_file`** is provided under the `grand` block, the Grand Canonical and Semi-Grand Canonical moves must be disabled to prevent stochastic fluctuations during the survey.
* **Single-Species Focus**: Only one solute (e.g., `Cr`) and one interstitial (e.g., `B`) should be defined. Parameters like `base_mu` and `c_target` are ignored in this mode as the process is not stochastic.

### 🚶 The Mapping Process
1. **Jump Logic**: The simulation performs a series of jumps within the provided void network. The `num_mc_steps` defines the total number of successful jumps the interstitial will attempt.
2. **Iterative Decoration**: At every site, the tool systematically replaces metal neighbors with the solute species. By default, it iterates from **0 to 6 nearest neighbors**.
3. **Barrier Calculation (`do_neb`)**: 
   * If `True`: The tool executes a **Climbing Image Nudged Elastic Band (CI-NEB)** calculation for every jump at every decoration level, recording both forward ($dE_f$) and reverse ($dE_r$) barriers.
   * If `False`: Only the total system energy is recorded at each site/decoration level, which is ideal for calculating **Segregation Energy ($\Delta E_{seg}$)** spectra without the kinetic overhead.


---

## 📏 Parameters & Suggestions

* **`jump_cutoff`**: Set to **2.15 Å** for this bulk Ni example (Octahedral-to-Octahedral). For custom grain boundaries, use the value suggested by the `void_generation` script.
* **`snapshot_every: 1`**: Set to 1 to ensure every unique decoration state and jump configuration is recorded for post-analysis.
* **`num_mc_steps`**: In this mode, the total simulation duration is calculated as $6 \times num\_mc\_steps$ to account for the iterative decoration steps.

---

## 💡 Usage Notes

### Determining the Spectrum
The resulting `spectral_log.csv` allows you to plot the interstitial energy as a function of the number of Cr neighbors. This "spectrum" reveals whether the solute acts as a trap or a repellant for the interstitial.

### Running the Example
```bash
python main.py --config config.yaml