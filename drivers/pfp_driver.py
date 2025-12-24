import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

def get_calculator():
    estimator = Estimator(model_version="v5.0.0", 
                          calc_mode=EstimatorCalcMode.PBE_PLUS_D3)
    return ASECalculator(estimator)