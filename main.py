import sys
import logging
from pathlib import Path

# Modular imports
from utils.config_parser import load_config
from utils.logger import setup_logging
from core.simulation import HybridSimulation, SpectralCollector

def main():
    # 1. Setup Logging
    logger = setup_logging()

    # 2. Load Configuration
    config_path = Path("config.yaml") 
    if not config_path.exists():
        logging.error(f"Configuration file {config_path} not found.")
        sys.exit(1)
        
    try:
        config = load_config(config_path)
        logging.info("Configuration loaded successfully.")
    except Exception as e:
        logging.error(f"Failed to parse configuration: {e}")
        sys.exit(1)
    
    # 3. Setup Workspace
    for folder in ["structures", "data", "data/neb"]:
        Path(folder).mkdir(parents=True, exist_ok=True)
    logging.info("Directory structure verified.")

    # 4. Initialize and Execute
    try:
        # Check the new dedicated spectral category
        spectral_config = config.get('spectral', {})
        is_spectral = spectral_config.get('enabled', False)

        if is_spectral:
            # Safety Check: Spectral mode requires grand canonical 'ingredients'
            # (voids_file and additives) even if random GC moves are disabled.
            if 'grand' not in config['ensembles']:
                logging.error("Spectral mode enabled, but 'ensembles: grand' section is missing.")
                sys.exit(1)
            
            logging.info("--- Mode: Spectral/NEB Mapping ---")
            # Inherits driver/calculator init from HybridSimulation
            sim = SpectralCollector(config) 
            sim.run_collection() 
        else:
            logging.info("--- Mode: Standard Hybrid MC/MD Simulation ---")
            sim = HybridSimulation(config)
            sim.run() 
        
        logging.info("Task completed successfully.")
        
    except KeyboardInterrupt:
        logging.warning("Simulation interrupted by user. Checkpoints saved.")
        sys.exit(0)
    except Exception as e:
        logging.exception(f"Fatal error during execution: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()