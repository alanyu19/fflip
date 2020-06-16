# -*- coding: utf-8 -*-


from fflip.chm.lipids.dppc import *
from fflip.chm.lipids.prpc import *
from fflip.chm.lipids.pc import *
from fflip.chm.lipids.pe import *
from fflip.chm.lipids.pg import *
from fflip.chm.lipids.pce import *


def parse_lipid(lipname):
    if lipname.lower() in ['dlpc', 'dppc', 'dmpc', 'popc', 'dopc', 'prpc', 'pc']:
        return pc
    if lipname.lower() in ['dppe', 'dope', 'dlpe', 'pope', 'pe']:
        return pe
    if lipname.lower()[-2:] == 'pg':
        return pg
    if lipname.lower()[-3:] == 'pce':
        return pce
