# -*- coding: utf-8 -*-

from fflip.drude.dlipid import *

# ******************************************************************************
# PC GROUPS:
drude_pc_charge_groups = list()
drude_pc_lj_groups = list()
# Please note that the CHARMM/DRUDE Force Field use Rmin/2 intead of sigma,
# and the unit of epsilon is kcal/mol


drude_pc_charge_groups.append(
    DrudeChargeGroup(
        id_=1,
        atom_groups=[
            ['C11'],
            ['C1']
        ],
        drude_particles=[
            ['DC11'],
            ['DC1']
        ],
        charges=[0.228, 0.026],
        alphas=[-1.642, -1.642],
        tholes=[0.862, 0.862],
        add_group=[True, False],
        add_alpha=[True, True],
        neighbors=[[1], []],
        atoms_same_charge=None,
        atoms_same_alpha=None
    )
)

drude_pc_charge_groups.append(
    DrudeChargeGroup(
        id_=2,
        atom_groups=[
            ['C2'], ['O21'], ['C21'], ['O22'], ['C22'],
            ['LPMA', 'LPMB'], ['LP1A'], ['LP1B']  # Virtual Sites (Lone Pairs)
        ],
        drude_particles=[
            ['DC2'], ['DO21'], ['DC21'], ['DO22'], ['DC22'],
            None, None, None
        ],
        charges=[
            0.202, 0.000, 0.697, 0.000, -0.206,
            -0.170, -0.349, -0.258
        ],
        alphas=[-1.797,  -0.732, -1.370, -0.904, -1.993, 0, 0, 0],
        tholes=[0.410, 0.601, 1.747, 0.565, 0.410, 0, 0, 0],
        add_group=[True, False, True, False, True, True, True, True],
        add_alpha=[True, True, True, True, True, False, False, False],
        neighbors=[[5], [], [4, 5, 6, 7], [], [2], [0, 2], [7, 2], [6, 2]],
        atoms_same_charge=[
            ['C3'], ['O31'], ['C31'], ['O32'], ['C32'],
            ['LPMC', 'LPMD'], ['LP1C'], ['LP1D']
        ],
        exclusion=[],
        atoms_same_alpha=[
            [('DC3', 'C3')], [('DO31', 'O31')], [('DC31', 'C31')],
            [('DO32', 'O32')], [('DC32', 'C32')], [], [], []
        ]
    )
)

drude_pc_charge_groups.append(
    DrudeChargeGroup(
        id_=3,
        atom_groups=[
            ['C3'], ['O31'], ['C31'], ['O32'], ['C32'],
            ['LPMC', 'LPMD'], ['LP1C'], ['LP1D']
        ],
        drude_particles=[
            ['DC3'], ['DO31'], ['DC31'], ['DO32'], ['DC32'],
            None, None, None
        ],
        charges=[
            0.086, 0.000, 0.697, 0.000, -0.206,
            -0.170, -0.349, -0.258
        ],
        alphas=[-1.797, -0.732, -1.370, -0.904, -1.993, 0, 0, 0],
        tholes=[0.410, 0.601, 1.747, 0.565, 0.410, 0, 0, 0],
        add_group=[False, False, False, False, False, False, False, False],
        add_alpha=[False, False, False, False, False, False, False, False],
        neighbors=[[], [], [], [], [], [], [], []],
        atoms_same_charge=None,
        atoms_same_alpha=None
    )
)

# drude_pc_lj_groups.append(
#     LJGroup(
#         id_=5,
#         atom_type_dict={
#             'ND3P2A': ['N'],
#             'CD33A': ['C13'],
#             'CD32A': ['C12'],
#             'PD1A': ['P'],
#             'OD2C2B': ['O13'],
#             'OD30B': ['O11'],
#             'CD32B': ['C1'],
#             'CD31C': ['C2'],
#             'OD30C': ['O21'],
#             'CD2O3A': ['C21'],
#             'OD2C3A': ['O22'],
#             'CD32C': ['C22']
#         }
#     )
# )

drude_pc_nbthole_groups = list()
drude_pc_nbthole_groups.append(
    NBTHOLEGroup(
        id_=4,
        nbthole_types={
            ('CD32B', 'OD2C2B'): 1.45
        }
    )
)
#  *****************************************************************************

# drude_pc = DrudeLipid(name="DPPC", charge_groups=drude_pc_charge_groups, lj_groups=drude_pc_lj_groups)
drude_dspc = DrudeLipid(
    name="DSPC", charge_groups=drude_pc_charge_groups,
    lj_groups=drude_pc_lj_groups, nbthole_groups=drude_pc_nbthole_groups
)

drude_dppc = DrudeLipid(
    name="DPPC", charge_groups=drude_pc_charge_groups,
    lj_groups=drude_pc_lj_groups, nbthole_groups=drude_pc_nbthole_groups
)

drude_dmpc = DrudeLipid(
    name="DMPC", charge_groups=drude_pc_charge_groups,
    lj_groups=drude_pc_lj_groups, nbthole_groups=drude_pc_nbthole_groups
)

drude_dlpc = DrudeLipid(
    name="DLPC", charge_groups=drude_pc_charge_groups,
    lj_groups=drude_pc_lj_groups, nbthole_groups=drude_pc_nbthole_groups
)

drude_popc = DrudeLipid(
    name="POPC", charge_groups=drude_pc_charge_groups,
    lj_groups=drude_pc_lj_groups, nbthole_groups=drude_pc_nbthole_groups
)
#  *****************************************************************************

match_lipid = {'dppc': drude_dppc, 'dmpc': drude_dmpc, 'dlpc': drude_dlpc, 'dspc': drude_dspc, 'popc': drude_popc}
