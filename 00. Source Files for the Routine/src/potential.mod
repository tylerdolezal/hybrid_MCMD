# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Define interatomic potentials
pair_style pfp_api v5.0.0 CRYSTAL_U0_PLUS_D3
pair_coeff * * species B Cr Ni
