# Example: SGC Equilibration of a Grain Boundary with Interstitials

This example demonstrates the final "Research-Production" workflow. It combines a **converged chemical potential ($\Delta\mu$)** with a **custom grain boundary geometry** to simulate the equilibrium segregation of solutes and the stability of interstitials.


## 🔬 Simulation Logic

At this stage, we transition from bulk calibration to localized interface physics. The simulation uses a pre-calculated $\Delta\mu$ to maintain a global chemical reservoir while allowing the grain boundary (GB) to naturally attract or repel species.

### 🛠️ Key Configuration Highlights
* **Converged $\Delta\mu$**: The `delta_mu` is set to **0.936**, a value previously calibrated to yield a $75/25$ Ni-Cr ratio at **300 K** for the **PFP** potential.
* **Custom Interface**: `use_custom_cell: True` loads a $\Sigma 5$ grain boundary (`POSCAR-custom`) that includes a pre-inserted Boron atom.
* **Interstitial Protection**: The randomization scheme is aware of the Boron's presence, ensuring the interstitial site is preserved while the metal lattice (Ni/Cr) is scrambled to the target composition.


## ⚙️ Simulation Setup

* **Ensembles**: Both **Canonical** and **Semi-Grand Canonical** are active. This allows for both the swapping of metal positions (Ordering) and the flipping of metal identities (Segregation).
* **Hybrid MD**: `hybrid_md` is enabled, running 20,000 steps of MD relaxation every 2,000 MC steps. This is vital for allowing the GB core to structurally reorganize as Chromium segregates to the interface.
* **Duration**: Set for 20,000 MC steps, providing a robust window for observing solute-interstitial interactions.

---

## 💡 Usage Notes

### Equilibrium Segregation
Because $\Delta\mu$ is fixed to the bulk-converged value, any deviation in Cr concentration at the grain boundary reflects the **Gibbsian Segregation** energy of the interface. If Cr prefers the GB, you will see a local concentration higher than 25% in that region.

### Interstitial Stability
With the Boron atom pre-inserted, you can monitor the local environment. Does Cr cluster around the Boron? Does the total energy of the system decrease as the Ni-Cr environment relaxes around the interstitial?


### Running the Example
```bash
python main.py --config config.yaml