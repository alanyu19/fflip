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
        add_alpha=[True, False, True, True, False],
        neighbors=[[2, 3], [4], [1], [4], []],
        atoms_same_charge=None,
        atoms_same_alpha=None
    )
)


drude_pc_charge_groups.append(
    DrudeChargeGroup(
        atom_groups=[
            ['C11', 'C1'],
            ['H11A', 'H11B', 'HA', 'HB'],
            ['P'],
            ['O13', 'O14'],
            ['O11', 'O12']
        ],
        drude_particles=[
            ['DC11', 'DC1'],
            None,
            ['DP'],
            ['DO13', 'DO14'],
            ['DO11', 'DO12']
        ],
        charges=[0.1985, 0.026, 1.191, -0.876, -0.470],
        alphas=[-1.642, 0, -0.974, -0.931, -0.901],
        tholes=[0.862, 2.098, 1.083, 0.181, 0.862],
        add_group=[True, True, True, True, True],
        add_alpha=[True, False, True, True, True],
        neighbors=[[4], [0], [3], [], [2]],
        atoms_same_charge=None,
        atoms_same_alpha=None
    )
)

drude_pc_charge_groups.append(
    DrudeChargeGroup(
        atom_groups=[
            ['C2'], ['HS'], ['O21'], ['C21'], ['O22'], ['C22'], ['H2R', 'H2S'],
            ['LP1A'], ['LP1B'], ['LPMA', 'LPMB']  # Virtual Sites (Lone Pairs)
        ],
        drude_particles=[
            ['DC2'], None, ['DO21'], ['DC21'], ['DO22'], ['DC22'], None,
            None, None, None
        ],
        charges=[
            0.202, 0.116, 0.000, 0.697, 0.000, -0.206, 0.069,
            -0.349, -0.258, -0.170
        ],
        alphas=[-1.797, 0, -0.732, -1.370, -0.904, -1.993, 0, 0, 0, 0],
        tholes=[0.410, 0, 0.601, 1.747, 0.565, 0.410, 0, 0, 0, 0],
        add_group=[
            True, False, False, True, False, True, True,
            True, True, True
        ],
        add_alpha=[
            True, False, True, True, True, True, False,
            False, False, False
        ],
        neighbors=[[1], [], [], [9, 5], [], [9], [5], [8, 3], [7, 3], [0]],
        # When add_charge_gtcnp is False but atoms_same_charge is not empty,
        # the later is used in neighbors/cooperators.
        # 'C3' is needed here (although the charge is not equal to C2)
        # since 'O31' will exchange charge with it
        # it will get skipped by the exclusion argument
        atoms_same_charge=[
            ['C3'], [], ['O31'], ['C31'], ['O32'], ['C32'], ['H2X', 'H2Y'],
            ['LP1C'], ['LP1D'], ['LPMC', 'LPMD']
        ],
        exclusion=[0],
        atoms_same_alpha=[
            [('DC3', 'C3')], [], [('DO31', 'O31')], [('DC31', 'C31')],
            [('DO32', 'O32')], [('DC32', 'C32')], [], [], [], []
        ]
    )
)

drude_pc_charge_groups.append(
    DrudeChargeGroup(
        atom_groups=[
            ['C3'], ['HX', 'HY'], ['O31'], ['C31'], ['O32'], ['C32'],
            ['H2X', 'H2Y'], ['LP1C'], ['LP1D'], ['LPMC', 'LPMD']
            # last three are Virtual Sites (Lone Pairs)
        ],
        drude_particles=[
            ['DC3'], None, ['DO31'], ['DC31'], ['DO32'], ['DC32'], None,
            None, None, None
        ],
        charges=[
            0.086, 0.116, 0.000, 0.697, 0.000, -0.206, 0.069,
            -0.349, -0.258, -0.170
        ],
        alphas=[-1.797, 0, -0.732, -1.370, -0.904, -1.993, 0, 0, 0, 0],
        tholes=[0.410, 0, 0.601, 1.747, 0.565, 0.410, 0, 0, 0, 0],
        add_group=[
            True, False, False, False, False, False, False,
            False, False, False],
        add_alpha=[
            False, False, False, False, False, False, False,
            False, False, False
        ],
        neighbors=[[1], [], [], [], [], [], [], [], [], []],
        atoms_same_charge=None,
        atoms_same_alpha=None
    )
)

#  *****************************************************************************

drude_pc = DrudeLipid(lipname="PC", charge_groups=drude_pc_charge_groups)

#  *****************************************************************************
