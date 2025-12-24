import logging
import sys

def setup_logging(level=logging.INFO):
    """
    Configures a dual-output logger that writes to both the console 
    and a 'simulation.log' file.
    """
    logger = logging.getLogger()
    
    # Avoid adding multiple handlers if setup_logging is called repeatedly
    if logger.hasHandlers():
        return logger

    logger.setLevel(level)

    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s', 
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # 1. Console Handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # 2. File Handler
    file_handler = logging.FileHandler('simulation.log', mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logging.info("--- Hybrid MC/MD Simulation Logger Initialized ---")
    return logger