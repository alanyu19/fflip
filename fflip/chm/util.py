# -*- coding: utf-8 -*-

import numpy as np
from fflip.chm.lipids.pc import *
from fflip.chm.lipids.pe import *
from fflip.chm.lipids.pg import *
from fflip.chm.lipids.pce import *
from fflip.chm.lipids.sm import *


def parse_lipid(lipname):
    if lipname.lower() in [
        'dlpc', 'dppc', 'dmpc', 'popc', 'dopc', 'prpc', 'pc'
    ]:
        return pc
    if lipname.lower() in ['dppe', 'dope', 'dlpe', 'pope', 'pe']:
        return pe
    if lipname.lower()[-2:] == 'pg':
        return pg
    if lipname.lower()[-3:] == 'pce':
        return pce
    if lipname.lower() in ['psm', 'ssm', 'sm']:
        return sm


def cosine_series(k, n, p, interval=5):
    """
    energy profile based on user provided:
    k: the force constant(s), can be array like or float
    n: the multiplicity(ies)
    p: the phase(s) in degrees
    return:
    angles, energies
    """
    energy = 0
    angles = np.arange(-180, 180, interval)
    angles = angles * np.pi / 180
    if type(k) == float or type(k) == int:
        k = [k]
        n = [n]
        p = [p]
    p = np.array(p)
    for k_, n_, p_ in zip(k, n, p):
        energy += k_ * (1 + np.cos(n_*angles - p_))
    return 180*angles/np.pi, energy
        
