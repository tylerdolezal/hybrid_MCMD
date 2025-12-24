import yaml
from pathlib import Path

def load_config(config_path="config.yaml"):
    """
    Loads and validates the simulation configuration from the finalized YAML structure.
    Returns the nested dictionary to maintain ensemble-awareness.
    """
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found at: {path}")

    with open(path, 'r') as f:
        try:
            config_data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML configuration: {e}")

    # --- 1. Validation of Section Presence ---
    # 'spectral' is now a top-level requirement alongside the others
    required_sections = ['simulation', 'system', 'ensembles', 'spectral']
    for section in required_sections:
        if section not in config_data:
            raise KeyError(f"Missing required configuration section: '{section}'")

    # --- 2. Deep Validation of Sub-keys ---
    check_map = {
        'simulation': ['num_mc_steps', 'temperature', 'potential_style'], # Cascade removed
        'system': ['use_custom_cell', 'composition'],
        'ensembles': ['canonical', 'semi_grand', 'grand'],
        'spectral': ['enabled', 'do_neb', 'jump_cutoff', 'solute'] # New section
    }

    for section, keys in check_map.items():
        for key in keys:
            if key not in config_data[section]:
                raise KeyError(f"Missing required key '{key}' in section '{section}'")

    # --- 3. Dependency Validation ---
    # Spectral mode requires the 'grand' section to be present for voids and additives
    if config_data['spectral']['enabled']:
        grand_cfg = config_data['ensembles']['grand']
        
        # We need the voids file to know the diffusion path
        if 'voids_file' not in grand_cfg:
            raise KeyError("Spectral mode is enabled but 'ensembles: grand: voids_file' is missing.")
            
        # We need additives to know what species is moving
        if 'additives' not in grand_cfg or not grand_cfg['additives']:
            raise KeyError("Spectral mode is enabled but 'ensembles: grand: additives' are missing.")

    return config_data