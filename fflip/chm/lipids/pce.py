# -*- coding: utf-8 -*-

from fflip.chm.lipid import *


# ******************************************************************************
# PC GROUPS:
pce_gs = list()
# Please note that the CHARMM force field use Rmin/2 intead of sigma,
# and the unit of epsilon is kcal/mol
pce_gs.append(
    CharmmGroup(
        num_atom_category=5,  # NTL, HL, CTL5, CTL2, HL
        atoms=[['N'],
               ['H13A', 'H13B', 'H13C', 'H14A', 'H14B', 'H14C', 'H15A', 'H15B',
                'H15C'],
               ['C13', 'C14', 'C15'],
               ['C12'],
               ['H12A', 'H12B']],
        charges=[-0.6, 0.25, -0.35, -0.1, 0.25],
        half_r_mins=[1.85, 0.7, 2.06, 2.10, 0.7],
        epsilons=[-0.20, -0.046, -0.0800, -0.0560, -0.046],
        add_charge_nbgroup=[True, True, True, True, False],
        atoms_same_charge=None,
        neighbors=[[2, 3], [4], [1], [4], []],
        cooperators=None,
        add_lj_nbgroup=[True, True, True, True, False],
        atoms_same_lj=[
            [], ['H12A', 'H12B'], [], ['C11', 'C1', 'C22', 'C3', 'C32'], []]
    )
)

pce_gs.append(
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
        add_charge_nbgroup=[True, True, True, True, True, True, False],
        atoms_same_charge=None,
        neighbors=[[2], [0], [3, 4], [4], [5], [6], []],
        cooperators=None,
        add_lj_nbgroup=[False, True, True, True, True, False, False],
        atoms_same_lj=[[], ['HA', 'HB', 'H2R', 'H2S', 'HX', 'HY', 'H2X', 'H2Y'],
                       [], [], [], [], []]
    )
)

pce_gs.append(
    CharmmGroup(
        num_atom_category=7,
        atoms=[
            ['C2'],  # CTL1 x
            ['HS'],  # HAL1 x
            ['O21'],  # OG301 x
            ['C21'],  # CTL2 x
            ['H1R', 'H1S'],  # HAL2 x
            ['C22'],  # CTL2 x
            ['H2R', 'H2S']  # HAL2 x
        ],
        charges=[0.13, 0.01, -0.56, 0.4, 0.01, -0.18, 0.09],
        half_r_mins=[2.275, 1.3200, 1.6500, 2.010, 1.3400, 2.010, 1.3400],
        epsilons=[-0.0200, -0.022, -0.1000, -0.0560, -0.028, -0.0560, -0.028],
        add_charge_nbgroup=[True, False, True, True, True, True, True],
        # When add_charge_nbgroup is False but atoms_same_charge is not empty,
        # the later is used in neighbors/cooperators.
        # 'C3' is needed here (although the charge is not equal to C2)
        # since 'O31' will exchange charge with it
        # it will get skipped by the exclusion argument
        atoms_same_charge=[
            ['C3'], [], ['O31'], ['C31'], [], ['C32'], ['H2X', 'H2Y']],
        exclusion=[0],
        neighbors=[[1], [], [0], [2, 5], [3], [3], [5]],
        cooperators=None,
        add_lj_nbgroup=[False, False, False, False, False, False, False],
        atoms_same_lj=None
    )
)

pce_gs.append(
    CharmmGroup(
        num_atom_category=7,
        atoms=[
            ['C3'],  # CTL2 x
            ['HX', 'HY'],  # HAL2 x
            ['O31'],  # OG301 x
            ['C31'],  # CTL2 x
            ['H1X', 'H1Y'],  # HAL2 x
            ['C32'],  # CTL2 x
            ['H2X', 'H2Y']  # HAL2 x
        ],
        charges=[0.08, 0.02, -0.56, 0.4, -0.02, -0.18, 0.09],
        half_r_mins=[2.010, 1.3400, 1.6500, 2.010, 1.3400, 2.010, 1.3400],
        epsilons=[-0.0560, -0.028, -0.1000, -0.0560, -0.028, -0.0560, -0.028],
        add_charge_nbgroup=[True, False, False, False, True, False, False],
        atoms_same_charge=None,
        neighbors=[[1], [], [], [], [3], [], []],
        cooperators=None,
        add_lj_nbgroup=[False, False, False, False, False, False, False],
        atoms_same_lj=None
    )
)

#  *****************************************************************************

pce = Lipid(charmm_group_list=pce_gs, name="PCE")

#  *****************************************************************************
