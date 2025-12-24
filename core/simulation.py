import numpy as np
import random
import importlib
import os
from ase.io import read, write
import core.thermo as thermo
import core.moves as moves
import utils.io as io
import logging
import pandas as pd

class HybridSimulation:
    def __init__(self, config):
        """
        Initializes simulation state, move tracking, and the physics driver.
        """
        self.config = config
        self.mc_step = 0
        self.accepted = 0
        self.rejected = 0
        self.temp = self.config['simulation']['temperature']
        self.is_spectral = config.get('spectral', {}).get('enabled', False)
        
        # 1. Dynamically load the potential driver (PFP, CHGNet, or EAM)
        self.calculator = self._init_driver()
        
        # 2. Initialize Chemical Potential Cache (Internalizes scaling logic)
        self.mu_cache = self._initialize_mu_cache()
        
        # Move tracking statistics
        self.move_stats = {
            'accepted': {m: 0 for m in ['swap', 'new_host', 'swap_ints', 'flip', 'insert', 'delete', 'shuffle', 'MD']},
            'rejected': {m: 0 for m in ['swap', 'new_host', 'swap_ints', 'flip', 'insert', 'delete', 'shuffle', 'MD']}
        }
        
        # Load system and energy history
        self.atoms, self.energies = self._setup_system()
        self.md_step_counter = 0  
        self.move_selector = moves.MoveSelector(config)

    def _initialize_mu_cache(self):
        """
        Calculates effective chemical potentials for Grand Canonical (GC) 
        and reads/sets delta_mu for Semi-Grand Canonical (SGC).
        """
        cache = {'additives': {}, 'delta_mu': 0.0}
        
        # Process Grand Canonical Additives
        # We check IF we are in spectral mode OR if grand is explicitly enabled
        grand_cfg = self.config['ensembles']['grand']
        if grand_cfg.get('enabled') or self.is_spectral:
            additives = grand_cfg.get('additives', {})
            for symbol, params in additives.items():
                # We still need the mu calculation for energy references
                eff_mu = thermo.get_effective_mu(
                    float(params['base_mu']), 
                    float(params.get('c_target')), 
                    self.temp
                )
                cache['additives'][symbol] = eff_mu
                logging.info(f"Species {symbol} registered for simulation (Spectral={self.is_spectral})")
        
        # Process Semi-Grand Canonical delta_mu
        if self.config['ensembles']['semi_grand']['enabled']:
            dm_input = self.config['ensembles']['semi_grand']['delta_mu']
            if isinstance(dm_input, str) and os.path.exists(dm_input):
                cache['delta_mu'] = float(np.loadtxt(dm_input))
            else:
                cache['delta_mu'] = float(dm_input)
                
        return cache

    def _init_driver(self):
        """Loads the specific potential driver."""
        style = self.config.get('potential_style', 'pfp').lower()
        try:
            driver_module = importlib.import_module(f"drivers.{style}_driver")
            return driver_module.get_calculator()
        except ImportError:
            raise ImportError(f"Driver drivers/{style}_driver.py not found.")

    def _setup_system(self):
        """Initializes system with the Modular Calculator."""

        if self.is_spectral:
            # FORCE standard behavior for Spectral Mode:
            # We ignore 'composition', 'randomize_initial', and 'crystal_shape'
            logging.info("Spectral mode active: Loading POSCAR-custom and bypassing randomization.")
            atoms = read("POSCAR-custom", format="vasp")
            return atoms, [0.0, 0.0]

        if self.config.get('continue_run'):
            atoms = read('POSCAR-1', format='vasp')
            energies = list(np.loadtxt('data/energies'))
            stats = io.parse_mc_statistics()
            self.mc_step = stats['steps_completed']
            self.move_stats['accepted'] = stats['accepted_counts']
            self.move_stats['rejected'] = stats['rejected_counts']
            return atoms, energies
        else:
            atoms, e0 = thermo.initialize_system(self.config, self.calculator)
            return atoms, [e0,e0]

    def run(self):
        """Main loop orchestrating Hybrid MD-MC steps."""
        total_steps = self.config['simulation']['num_mc_steps']
        md_config = self.config['simulation']['hybrid_md']

        for step in range(self.mc_step + 1, total_steps + 1):
            self.mc_step = step
            
            # Check MD toggle logic
            if md_config['enabled'] and self.md_step_counter >= md_config['interval']:
                self.execute_md_step()
                self.md_step_counter = 0
            else:
                self.execute_mc_step()
                self.md_step_counter += 1

            if step % self.config['simulation'].get('snapshot_every', 200) == 0:
                io.save_snapshot(self)

    def execute_mc_step(self):
        """Evaluates a Metropolis move using cached mu values."""
        weights = self.move_selector.get_adaptive_weights(self.atoms, self.move_stats)
        move_list = ['flip', 'insert', 'delete', 'swap', 'new_host', 'shuffle', 'swap_ints']
        move_type = random.choices(move_list, weights=weights, k=1)[0]
        
        swap_pairs = self.move_selector.select_atoms(self.atoms, move_type, self.move_stats)
        
        # Pass mu_cache to thermo to avoid recalculating scaled mu
        prob, new_atoms, new_energy, new_voids = thermo.calculate_energy_change(
            self.atoms, self.energies[-1], swap_pairs, move_type, 
            self.config, self.calculator, self.mu_cache
        )
        
        if random.random() < prob:
            self.atoms = new_atoms
            self.energies.append(new_energy)
            self.move_stats['accepted'][move_type] += 1
            self.accepted += 1
            if new_voids is not None:
                write(self.config['ensembles']['grand']['voids_file'], new_voids, format='vasp', sort=True)
        else:
            self.energies.append(self.energies[-1])
            self.move_stats['rejected'][move_type] += 1
            self.rejected += 1

    def execute_md_step(self):
        """Dispatches structure to the MD routine."""
        print(f"--- Step {self.mc_step}: Running MD ---")
        md_config = self.config['simulation']['hybrid_md']
        
        # Standardized MD call using global temperature
        new_atoms, new_energy = thermo.run_md(
            self.atoms, 
            self.calculator, 
            self.temp, 
            md_config['steps']
        )
        
        self.atoms = new_atoms
        self.energies.append(new_energy)
        self.move_stats['accepted']['MD'] += 1
        self.accepted += 1

class SpectralCollector(HybridSimulation):
    def __init__(self, config):
        """
        Initializes the collector by inheriting all physics drivers 
        and setup from HybridSimulation.
        """
        super().__init__(config)
        self.solute = self.config['spectral'].get('solute', 'Cr')
        self.run_neb = self.config['spectral'].get('do_neb', False)
        self.log_path = "data/spectral_log.csv"

        # Initialize spectral log
        self._init_spectral_log()

    def run_collection(self):
        logging.info("Starting Random-Walk Diffusion Scan...")
        
        # 1. Setup
        all_void_atoms = read(self.config['ensembles']['grand']['voids_file'])
        symbol = list(self.mu_cache['additives'].keys())[0]
        
        # 2. Initial Reference Site
        # Place at the first void and relax
        self.atoms += thermo.Atoms(symbol, positions=[all_void_atoms[0].position])
        int_idx = len(self.atoms) - 1
        self.atoms, _ = thermo.relax_config(self.atoms, 'insert', self.calculator)
        
        # write the poscar
        write("POSCAR-0", self.atoms, format='vasp', direct=True, sort=True)

        last_void_idx = 0 # Track where we just came from

        for step in range(self.config['simulation']['num_mc_steps']):
            logging.info(f"--- Diffusion Jump {step+1} ---")
            
            current_pos = self.atoms[int_idx].position
            jump_cutoff = self.config['spectral'].get('jump_cutoff', 3.0)
            
            # 3. Identify all valid neighbors (excluding the exact site we are on)
            valid_targets = []
            for i, v in enumerate(all_void_atoms):
                dist = np.linalg.norm(v.position - current_pos)
                # Filter: within cutoff, not current site, and ideally not the immediate last site
                if 0.1 < dist <= jump_cutoff and i != last_void_idx:
                    valid_targets.append(i)

            # Fallback: if trapped, allow jumping back to the last site
            if not valid_targets:
                valid_targets = [i for i, v in enumerate(all_void_atoms) 
                                 if 0.1 < np.linalg.norm(v.position - current_pos) <= jump_cutoff]

            if not valid_targets:
                logging.warning("Atom is trapped with no neighbor voids. Ending scan.")
                break

            # 4. Randomly pick a neighbor to prevent "teleporting" through the list
            target_idx = random.choice(valid_targets)
            target_void = all_void_atoms[target_idx]
            
            # Execute Jump
            initial_pos = current_pos.copy()
            self.atoms[int_idx].position = target_void.position
            
            # 5. Decoration Sweep & NEB
            for n in range(7):
                decorated_system = self.atoms.copy()
                if n > 0:
                    decorated_system = thermo.decorate_environment(
                        decorated_system, self.solute, n, int_idx, self.config
                    )
                
                final_system, final_energy = thermo.relax_config(
                    decorated_system, 'spectral', self.calculator
                )
                
                dEf, dEr = 0.0, 0.0
                if self.run_neb:
                    neb_start = final_system.copy()
                    neb_start[int_idx].position = initial_pos
                    dEf, dEr = thermo.run_neb_calculation(neb_start, final_system, self._init_driver, symbol)

                io.record_spectral_entry(final_energy, dEf, dEr, final_system, n, self.solute)

                # write the poscar
                write(f"POSCAR-1", final_system, format='vasp', direct=True, sort=True)

            # Update tracking for next jump
            last_void_idx = target_idx

    def _init_spectral_log(self):
        """Initializes the CSV headers for spectral mapping."""
        columns = ["energy (eV)", "dE_forward (eV)", "dE_reverse (eV)", "x (Å)", "y (Å)", "z (Å)", f"n_{self.solute}"]
        df = pd.DataFrame(columns=columns)
        df.to_csv(self.log_path, index=False)