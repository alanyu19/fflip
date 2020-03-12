# -*- coding: utf-8 -*-

from fflip.drude.dlipid import *


# ******************************************************************************
# PC GROUPS:
drude_pc_charge_groups = list()
# Please note that the CHARMM/DRUDE Force Field use Rmin/2 intead of sigma,
# and the unit of epsilon is kcal/mol

drude_pc_charge_groups.append(
    DrudeChargeGroup(
        atom_groups=[
            ['N'],
            ['H13A', 'H13B', 'H13C', 'H14A', 'H14B',
             'H14C', 'H15A', 'H15B', 'H15C'],
            ['C13', 'C14', 'C15'],
            ['C12'],
            ['H12A', 'H12B']
        ],
        drude_particles=[
            ['DN'],
            None,
            ['DC13', 'DC14', 'DC15'],
            ['DC12'],
            None
        ],
        charges=[0.688, 0.122, -0.288, -0.166, 0.122],
        alphas=[-0.829, 0, -1.793, -1.393, 0],
        tholes=[0.793, 0, 1.099, 0.862, 0],
        add_group=[True, True, True, True, False],
        neighbors=[[2, 3], [4], [1], [4], []],
        atoms_same_charge=None,
        atoms_same_alpha=None
    )
)

"""
pc_gs.append(
    CharmmGroup(
        num_atom_category=7,  # CTL2, HAL2, PL, 02L, OSLP, CTL2, HAL2
        atoms=[['C11'],
               ['H11A', 'H11B'],
               ['P'],
               ['O13', 'O14'],
               ['O11', 'O12'],
               ['C1'],
               ['HA', 'HB']],
        charges=[-0.08, 0.09, 1.50, -0.78, -0.57, -0.08, 0.09],
        half_r_mins=[2.010, 1.3400, 2.15, 1.70, 1.6500, 2.010, 1.3400],
        epsilons=[-0.0560, -0.028, -0.585, -0.12, -0.1000, -0.0560, -0.028],
        add_charge_gtcnp=[True, True, True, True, True, True, False],
        atoms_same_charge=None,
        neighbors=[[2], [0], [3, 4], [4], [5], [6], []],
        cooperators=None,
        add_lj_gtcnp=[False, True, True, True, True, False, False],
        atoms_same_lj=[[], ['HA', 'HB', 'H2R', 'H2S', 'HX', 'HY', 'H2X', 'H2Y'],
                       [], [], [], [], []]
    )
)

pc_gs.append(
    CharmmGroup(
        num_atom_category=7,  # CTL1, HAL1, OSL, CL, OBL, CTL2, HAL2
        atoms=[
            ['C2'], ['HS'], ['O21'], ['C21'], ['O22'], ['C22'], ['H2R', 'H2S']
        ],
        charges=[0.17, 0.09, -0.49, 0.90, -0.63, -0.22, 0.09],
        half_r_mins=[2.275, 1.3200, 1.6500, 2.00, 1.70, 2.010, 1.3400],
        epsilons=[-0.0200, -0.022, -0.1000, -0.0700, -0.12, -0.0560, -0.028],
        add_charge_gtcnp=[True, False, True, True, True, True, True],
        # When add_charge_gtcnp is False but atoms_same_charge is not empty,
        # the later is used in neighbors/cooperators.
        # 'C3' is needed here (although the charge is not equal to C2)
        # since 'O31' will exchange charge with it
        # it will get skipped by the exclusion argument
        atoms_same_charge=[
            ['C3'], [], ['O31'], ['C31'], ['O32'], ['C32'], ['H2X', 'H2Y']],
        exclusion=[0],
        neighbors=[[1], [], [0], [2, 4, 5], [3], [3], [5]],
        cooperators=None,
        add_lj_gtcnp=[True, True, True, True, True, False, False],
        atoms_same_lj=[[], [], ['O31'], ['C31'], ['O32'], [], []]
    )
)

pc_gs.append(
    CharmmGroup(
        num_atom_category=7,  # CTL2, HAL2, OSL, CL, OBL, CTL2, HAL2
        atoms=[['C3'], ['HX', 'HY'], ['O31'], ['C31'], ['O32'],
               ['C32'], ['H2X', 'H2Y']],
        charges=[0.08, 0.09, -0.49, 0.90, -0.63, -0.22, 0.09],
        half_r_mins=[2.275, 1.3200, 1.6500, 2.00, 1.70, 2.010, 1.3400],
        epsilons=[-0.0200, -0.022, -0.1000, -0.0700, -0.12, -0.0560, -0.028],
        add_charge_gtcnp=[True, False, False, False, False, False, False],
        atoms_same_charge=None,
        neighbors=[[1], [], [], [], [], [], []],
        cooperators=None,
        add_lj_gtcnp=[False, False, False, False, False, False, False],
        atoms_same_lj=None
    )
)
"""

#  *****************************************************************************

drude_pc = DrudeLipid(lipname="PC", charge_groups=drude_pc_charge_groups)

#  *****************************************************************************
