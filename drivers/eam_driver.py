from ase.calculators.lammpslib import LAMMPSlib

def get_calculator():
    # Example for Ni-Cr EAM potential
    cmds = ["pair_style eam/alloy",
            "pair_coeff * * NiCr.eam.alloy Ni Cr"]
    return LAMMPSlib(lmpcmds=cmds, atom_types={'Ni': 1, 'Cr': 2})