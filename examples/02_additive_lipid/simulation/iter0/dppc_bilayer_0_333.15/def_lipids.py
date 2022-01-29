# -*- coding: utf-8 -*-

from fflip.chm.lipid import *


# ******************************************************************************
# PC GROUPS:
pc_gs = list()
# Please note that the CHARMM force field use Rmin/2 intead of sigma,
# and the unit of epsilon is kcal/mol
pc_gs.append(
    CharmmGroup(
        num_atom_category=7,  # CTL1, HAL1, OSL, CL, OBL, CTL2, HAL2
        atoms=[
            ['C2'], ['HS'], ['O21'], ['C21'], ['O22'], ['C22'], ['H2R', 'H2S']
        ],
        charges=[0.17, 0.09, -0.49, 0.90, -0.63, -0.22, 0.09],
        half_r_mins=[2.275, 1.3200, 1.6500, 2.00, 1.70, 2.010, 1.3400],
        epsilons=[-0.0200, -0.022, -0.1000, -0.0700, -0.12, -0.0560, -0.028],
        add_charge_group=[True, False, True, True, True, True, True],
        # When add_charge_nbgroup is False but atoms_same_charge is not empty,
        # the later is used in neighbors/cooperators.
        # 'C3' is needed here (although the charge is not equal to C2)
        # since 'O31' will exchange charge with it
        # it will get skipped by the exclusion argument
        atoms_same_charge=[
            ['C3'], [], ['O31'], ['C31'], ['O32'], ['C32'], ['H2X', 'H2Y']],
        neighbors=[[1], [], [0], [2, 4, 5], [3], [3], [5]],
        add_lj_group=[False, False, True, False, True, False, False],
        atoms_same_lj=[[], [], ['O31'], ['C31'], ['O32'], [], []],
        exclusion=[0]
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
        add_charge_group=[True, False, False, False, False, False, False],
        atoms_same_charge=None,
        neighbors=[[1], [], [], [], [], [], []],
        add_lj_group=[False, False, False, False, False, False, False],
        atoms_same_lj=None
    )
)

#  *****************************************************************************

# name here is different from lipname in targets or sim/potential templates.
# this name is for a category of lipids.
pc_demo = Lipid(charmm_group_list=pc_gs, name="PC")

#  *****************************************************************************

# you also need to add 'dmpc' as pc_demo in match_lipid
# if you have one system containing 'dmpc'
match_lipid = {'dppc': pc_demo, 'dspc': pc_demo, 'dlpc': pc_demo, 'dmpc': pc_demo,
               'popc': pc_demo, 'sopc': pc_demo}
