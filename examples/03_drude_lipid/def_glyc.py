# -*- coding: utf-8 -*-

import sys
import nlopt
from scipy.optimize import minimize

from fflip.drude.dlipid import *
from fflip.omm.util import get_md_options as gmd
from fflip.ljpme.modelcomp import *
import fflip.omm.energiesforces as ef
# from rflow.observables import *
# from rflow.trajectory import *
# from fflip.chm import *
# from fflip.drude import *
# import rflow.observables as ro
from matplotlib import pyplot as plt


# ******************************************************************************
dcg1 = DrudeChargeGroup(
    id_=1,
    atom_groups=[
        ['CP1'],
        # ['HP11', 'HP12', 'HP13'],
        # ['P1'],
        # ['OP3', 'OP4'],
        # ['OP1', 'OP2'],
        ['CP2'],
        # ['H31', 'H32']
    ],
    drude_particles=[
        ['DCP1'],
        # None,
        # ['DP1'],
        # ['DOP3', 'DOP4'],
        # ['DOP1', 'DOP2'],
        ['DCP2'],
        # None
    ],
    # charges=[0.202, 0.026, 1.192, -0.856, -0.520, 0.228, 0.026],
    # alphas=[-1.642, 0, -0.974, -0.931, -0.901, -1.642, 0],
    # tholes=[0.862, 0, 2.098, 1.083, 0.181, 0.862, 0],
    # add_group=[True, True, True, True, True, True, True],
    # add_alpha=[True, False, True, True, True, True, False],
    # neighbors=[[2, 3, 4], [0], [3, 4], [2], [2], [2, 3, 4], [5]],
    charges=[0.202, 0.228],
    alphas=[-1.642, -1.642],
    tholes=[0.862, 0.862],
    add_group=[True, False],
    add_alpha=[True, True],
    neighbors=[[1], []],
    atoms_same_charge=None,
    atoms_same_alpha=None
)

dcg2 = DrudeChargeGroup(
    id_=2,
    atom_groups=[
        ['C2'],
        # ['H21'],
        ['O8'],
        ['C9'],
        ['O10'],
        ['C11'],
        # ['H111', 'H112', 'H113'],
        ['LP8A', 'LP8B'],
        ['LP1A'],
        ['LP1B'],
    ],
    drude_particles=[
        ['DC2'],
        # None,
        ['DO8'],
        ['DC9'],
        ['DO10'],
        ['DC11'],
        # None,
        None,
        None,
        None
    ],
    # charges=[0.202, 0.116, 0.000, 0.697, 0.000, -0.275, 0.069, -0.170, -0.349, -0.258],
    # alphas=[-1.797, 0, -0.732, -1.370, -0.904, -1.993, 0, 0, 0, 0],
    # tholes=[0.410, 0, 0.601, 1.747, 0.565, 0.410, 0, 0, 0, 0],
    # add_group=[True, False, False, True, False, True, True, True, True, True],
    # add_alpha=[True, False, True, True, True, True, False, False, False, False],
    # neighbors=[[1], [], [], [9, 5], [], [9], [5], [8, 3], [7, 3], [0]],
    # atoms_same_charge=[
    #     ['C1'], [], ['O4'], ['C5'], ['O6'], ['C7'], ['H71', 'H72', 'H73'],
    #     ['LP4A', 'LP4B'], ['LP6A'], ['LP6B']
    # ],
    # exclusion=[0],
    # atoms_same_alpha=[
    #     [('DC1', 'C1')], [], [('DO4', 'O4')], [('DC5', 'C5')],
    #     [('DO6', 'O6')], [('DC7', 'C7')], [], [], [], []
    # ]
    charges=[0.202, 0.000, 0.697, 0.000, -0.275, -0.170, -0.349, -0.258],
    alphas=[-1.797, -0.732, -1.370, -0.904, -1.993, 0, 0, 0],
    tholes=[0.410, 0.601, 1.747, 0.565, 0.410, 0, 0, 0],
    add_group=[True, False, True, False, True, True, True, True],
    add_alpha=[True, True, True, True, True, False, False, False],
    neighbors=[[5], [], [4, 5, 6, 7], [], [2], [0, 2], [7, 2], [6, 2]],
    atoms_same_charge=[
        ['C1'], ['O4'], ['C5'], ['O6'], ['C7'],
        ['LP4A', 'LP4B'], ['LP6A'], ['LP6B']
    ],
    exclusion=[],
    atoms_same_alpha=[
        [('DC1', 'C1')], [('DO4', 'O4')], [('DC5', 'C5')],
        [('DO6', 'O6')], [('DC7', 'C7')], [], [], []
    ]
)

dcg3 = DrudeChargeGroup(
    id_=3,
    atom_groups=[
        ['C1'],
        # ['H11', 'H12'],
        ['O4'],
        ['C5'],
        ['O6'],
        ['C7'],
        # ['H71', 'H72', 'H73'],
        ['LP4A', 'LP4B'],
        ['LP6A'],
        ['LP6B'],
    ],
    drude_particles=[
        ['DC1'],
        # None,
        ['DO4'],
        ['DC5'],
        ['DO6'],
        ['DC7'],
        # None,
        None,
        None,
        None
    ],
    charges=[0.202, 0.000, 0.697, 0.000, -0.275, -0.170, -0.349, -0.258],
    alphas=[-1.797, -0.732, -1.370, -0.904, -1.993, 0, 0, 0],
    tholes=[0.410, 0.601, 1.747, 0.565, 0.410, 0, 0, 0],
    add_group=[False, False, False, False, False, False, False, False],
    add_alpha=[False, False, False, False, False, False, False, False],
    neighbors=[[], [], [], [], [], [], [], []],
    atoms_same_charge=None,
    atoms_same_alpha=None
)

glyp_charge_groups = list()
glyp_lj_groups = list()
glyp_nbthole_groups = list()
# Please note that the CHARMM/DRUDE Force Field use Rmin/2 intead of sigma,
# and the unit of epsilon is kcal/mol

glyp_charge_groups += [dcg1, dcg2, dcg3]

#glyp_nbthole_groups.append(
#    NBTHOLEGroup(
#        id_=4,
#        nbthole_types={
#            ('CD32B', 'OD2C2B'): 1.45
#        }
#    )
#)

drude_glyc = DrudeLipid(
     name="GLYC", charge_groups=glyp_charge_groups,
     lj_groups=glyp_lj_groups, nbthole_groups=glyp_nbthole_groups
)

