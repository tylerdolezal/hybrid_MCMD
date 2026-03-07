import random
import numpy as np
from ase.neighborlist import NeighborList, natural_cutoffs

class MoveSelector:
    def __init__(self, config):
        """
        Initializes the MoveSelector with configuration parameters.
        """
        self.config = config
        self.interstitials = ["B", "C", "Cl", "O", "N", "H"]
        
        # Access additives from the new nested Grand Canonical path
        grand_config = self.config['ensembles']['grand']
        if grand_config['enabled'] and 'additives' in grand_config:
            add_species = list(grand_config['additives'].keys())
            self.interstitials = list(set(self.interstitials) | set(add_species))

    def get_adaptive_weights(self, atoms, stats):
        """
        Allocates equal weights to all moves that are both enabled 
        in the configuration and valid for the current system state.
        """
        try:
            symbols = atoms.get_chemical_symbols()
        except Exception:
            symbols = []

        # 1. System state analysis
        N_B = sum(symbols.count(sym) for sym in self.interstitials)
        metal_symbols = [sym for sym in symbols if sym not in self.interstitials]
        diverse_metal = len(set(metal_symbols)) > 1
        diverse_ints = len(set([sym for sym in symbols if sym in self.interstitials])) > 1
        
        ens_config = self.config['ensembles']
        active_moves = []

        # 2. Check validity for each move type
        # Semi-grand Ensemble
        if ens_config['semi_grand'].get('enabled', False):
            active_moves.append('flip')

        # Grand Canonical Ensemble
        if ens_config['grand'].get('enabled', False):
            active_moves.append('insert')
            if N_B > 0:
                active_moves.append('delete')

        # Canonical Ensemble
        if ens_config['canonical'].get('enabled', False):
            if diverse_metal:
                active_moves.append('swap')
            if N_B > 0:
                active_moves.extend(['new_host', 'shuffle'])
                if diverse_ints:
                    active_moves.append('swap_ints')

        # 3. Calculate equal weights
        weights = {m: 0.0 for m in ['flip', 'insert', 'delete', 'swap', 'new_host', 'shuffle', 'swap_ints']}
        
        if active_moves:
            equal_weight = 1.0 / len(active_moves)
            for move in active_moves:
                weights[move] = equal_weight

        # 4. Return ordered list
        return [
            weights['flip'], weights['insert'], weights['delete'],
            weights['swap'], weights['new_host'], weights['shuffle'],
            weights['swap_ints']
        ]

    def select_atoms(self, system, move_type, stats):
        """Dispatcher for selecting atom indices."""
        mc_step = stats.get('mc_step', 0)
        
        if move_type == 'swap':
            return self._select_swap_pairs(system)
        elif move_type == 'swap_ints':
            return self._select_swap_int_pairs(system)
        elif move_type == 'flip':
            return self._select_flip_atom(system)
        elif move_type == 'new_host':
            return self._select_host_diffusion(system, mc_step)
        elif move_type == 'shuffle':
            return self._select_shuffle_pairs(system)
        return None

    def _get_neighbors(self, system, idx, multiplier=1.0):
        """Finds valid neighbors using natural_cutoffs (no freeze threshold)."""
        cutoffs = natural_cutoffs(system, multiplier=multiplier)
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(system)
        indices, _ = nl.get_neighbors(idx)
        return list(indices)

    def _select_swap_pairs(self, system):
        """Selects two metal atoms of different species to exchange."""
        metal_indices = [i for i, a in enumerate(system) if a.symbol not in self.interstitials]
        if not metal_indices: return None
        
        host = random.choice(metal_indices)
        
        others = [i for i in metal_indices if system[i].symbol != system[host].symbol]
        if not others: return None
        return [(host, random.choice(others))]

    def _select_swap_int_pairs(self, system):
        """Selects two interstitial atoms to exchange."""
        int_indices = [i for i, a in enumerate(system) if a.symbol in self.interstitials]
        if not int_indices: return None

        # Select two distinct interstitials
        int1 = random.choice(int_indices)
        others = [i for i in int_indices if i != int1]
        if not others: return None
        int2 = random.choice(others)
        return [(int1, int2)]

    def _select_host_diffusion(self, system, mc_step):
        """
        Implements long-range diffusion by selecting metal hosts 
        outside the immediate 1NN shell.
        """
        int_indices = [i for i, a in enumerate(system) if a.symbol in self.interstitials]
        if not int_indices: 
            return None
        
        # Pick the interstitial to move
        target_int = random.choice(int_indices)
        
        # 1. Identify the 1NN "Forbidden" Shell
        # We use your existing neighbor logic to find who is too close
        forbidden_neighbors = set(self._get_neighbors(system, target_int))
        # Also include the current 'host' index if it's stored in the system
        forbidden_neighbors.add(target_int) 

        # 2. Identify all valid Metal Hosts
        # We define 'hosts' as any atom that is NOT an interstitial
        all_metals = [i for i, a in enumerate(system) if a.symbol not in self.interstitials]
        
        # 3. Filter: Any metal host NOT in the forbidden shell
        potential_hosts = [i for i in all_metals if i not in forbidden_neighbors]

        if not potential_hosts: 
            return None
            
        # Pick a random distant host
        target_host_idx = random.choice(potential_hosts)
        
        return [(target_int, target_host_idx)]

    def _select_flip_atom(self, system):
        """Decides between Local vs Global flips."""
        int_indices = [i for i, a in enumerate(system) if a.symbol in self.interstitials]
        choices = self.config['ensembles']['semi_grand'].get('metal_library', ['Ni', 'Cr'])
        
        if int_indices and random.random() < 0.5: # Local flip
            candidates = set()
            for idx in int_indices:
                for n in self._get_neighbors(system, idx):
                    if system[n].symbol in choices: candidates.add(n)
            candidate_indices = list(candidates)
        else: # Global flip
            candidate_indices = [i for i, a in enumerate(system) if a.symbol in choices]

        if not candidate_indices: return None
        return [(random.choice(candidate_indices), None)]

    def _select_shuffle_pairs(self, system):
        """Simplified global shuffling around interstitials."""
        int_indices = [i for i, a in enumerate(system) if a.symbol in self.interstitials]
        if not int_indices: return None
        
        target_int = random.choice(int_indices)
        neighbors = self._get_neighbors(system, target_int)
        if not neighbors: return None
        
        neighbor_idx = random.choice(neighbors)
        others = [i for i in range(len(system)) if system[i].symbol not in self.interstitials 
                  and system[i].symbol != system[neighbor_idx].symbol]
        
        if not others: return None
        return [(neighbor_idx, random.choice(others))]