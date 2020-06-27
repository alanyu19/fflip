# -*- coding: utf-8 -*-

from fflip.chm.lipid import *


# ******************************************************************************
# Sphingomyelin GROUPS:
sm_gs = list()
# Please note that the CHARMM force field use Rmin/2 intead of sigma,
# and the unit of epsilon is kcal/mol
sm_gs.append(
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
        add_charge_gtcnp=[True, True, True, True, False],
        atoms_same_charge=None,
        neighbors=[[2, 3], [4], [1], [4], []],
        cooperators=None,
        add_lj_gtcnp=[True, True, True, True, False],
        atoms_same_lj=[
            [], ['H12A', 'H12B'], [], ['C11', 'C1', 'C22', 'C3', 'C32'], []]
    )
)

sm_gs.append(
    CharmmGroup(
        num_atom_category=7,  # CTL2, HAL2, PL, 02L, OSLP, CTL2, HAL2
        atoms=[['C11'],
               ['H11A', 'H11B']
               ['P'],
               ['O13', 'O14'],
               ['O11', 'O12'],
               ['C1S'],
               ['H1S', 'H1T']],
        charges=[-0.08, 0.09, 1.50, -0.78, -0.57, -0.08, 0.09],
        half_r_mins=[2.010, 1.3400, 2.15, 1.70, 1.6500, 2.010, 1.3400],
        epsilons=[-0.0560, -0.028, -0.585, -0.12, -0.1000, -0.0560, -0.028],
        add_charge_gtcnp=[True, True, True, True, True, True, False],
        atoms_same_charge=None,
        neighbors=[[2], [0], [3, 4], [4], [5], [6], []],
        cooperators=None,
        add_lj_gtcnp=[False, True, True, True, True, False, False],
        # TODO: GO BACK TO CHANGE THIS PART
        atoms_same_lj=[[], ['HA', 'HB', 'H2R', 'H2S', 'HX', 'HY', 'H2X', 'H2Y'],
                       [], [], [], [], []]
    )
)

sm_gs.append(
    CharmmGroup(
        num_atom_category=7,  # CTL1, HAL1, NHL, H, C, O, HAL2
        atoms=[
            ['C2S'], ['H2S'], ['NF'], ['HNF'],
            ['C1F'], ['OF'], ['C2F'], ['H2F', 'H2G']
        ],
        charges=[0.30, 0.05, -0.70, 0.35, 0.55, -0.60, -0.07, 0.06],
        half_r_mins=[
            2.275, 1.3200, 1.8500, 0.2245,
            2.0000, 1.7000, 2.010, 1.3400],
        epsilons=[-0.0200, -0.0220, -0.2000, -0.0460,
                  -0.1100, -0.1200, -0.0560, -0.0280],
        add_charge_gtcnp=[True, True, True, True, True, True, True, True],
        # When add_charge_gtcnp is False but atoms_same_charge is not empty,
        # the later is used in neighbors/cooperators.
        # 'C3' is needed here (although the charge is not equal to C2)
        # since 'O31' will exchange charge with it
        # it will get skipped by the exclusion argument
        atoms_same_charge=None,
        neighbors=[[2], [0], [0, 4], [2], [2, 6], [4], [4], [6]],
        cooperators=None,
        add_lj_gtcnp=[False, False, False, False, False, False, False, False],
        atoms_same_lj=None  # since we don't change LJ, we don't need this
    )
)

#  *****************************************************************************

sm = Lipid(charmm_group_list=sm_gs, lipname="SM")

#  *****************************************************************************
