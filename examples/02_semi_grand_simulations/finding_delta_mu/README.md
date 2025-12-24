# Example: Semi-Grand Canonical (SGC) Calibration

This configuration is specifically designed to determine the **Target Chemical Potential Difference ($\Delta\mu$)** required to achieve a specific solute concentration in the bulk lattice. This calibration serves as the foundation for defining SGC reservoirs in more complex simulations.



## 🔬 Calibration Logic

The primary goal of this setup is to map a target mole fraction ($x_{solute}$) to a specific $\Delta\mu$. Note that the chemical potential is a function of the target chemistry and thermal state: **$\Delta\mu = f(x_{solute}, T)$**.

### 🛠️ Key Configuration Adjustments
* **`hybrid_md: enabled: False`**: Lattice relaxation is disabled to isolate the energetic cost of chemical species identity flips.
* **`num_mc_steps: 10000`**: The step count is reduced compared to full equilibration runs, as chemical identity swaps often converge faster than mechanical relaxations.
* **`snapshot_every: 100`**: Frequency is increased to closely monitor the trajectory of the composition changes.
* **`delta_mu: -0.25`**: A float value is provided as an initial guess based on the **transmutation energy** of a single host atom (e.g., Ni) into a solute atom (e.g., Cr).

## ⚙️ Simulation Setup

* **Ensemble**: **Semi-Grand Canonical** is enabled to allow metal identity flips (e.g., $Ni \leftrightarrow Cr$).
* **System**: Generates a standard $8 \times 8 \times 8$ FCC supercell.
* **Initial Scramble**: `randomize_initial: True` ensures the simulation starts from a disordered state, testing the ability of $\Delta\mu$ to drive the system to the target concentration.
* **Potential**: Utilizes the **PFP** machine-learning potential for high-fidelity transmutation energies.

---

## 💡 Usage Notes

### How to Calibrate
1. Set the target composition in `composition: {Ni: 0.85, Cr: 0.15}`.
2. Define the temperature (e.g., `temperature: 300`).
3. Run the simulation with your initial `delta_mu` guess.
4. Analyze the `simulation.log` to see the final average concentration of Cr.
5. If the final concentration is too low, decrease $\Delta\mu$ (make it more negative for Ni $\rightarrow$ Cr flips); if too high, increase it.

### ⚠️ Be Advised
* **Functional Dependency**: Because $\Delta\mu$ is dependent on the state point ($x_{solute}, T$), **you must re-calibrate $\Delta\mu$ if you change either the target composition or the simulation temperature**. 
* **Potential Specificity**: Each interatomic potential (PFP, CHGNet, EAM) interprets the "cost" of an atom flip differently. **Each potential will converge onto its own $\Delta\mu$**. You cannot reuse a $\Delta\mu$ value found using EAM for a simulation driven by PFP.

### Running the Example
```bash
python main.py --config config.yaml