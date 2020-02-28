# -*- coding: utf-8 -*-


from fflip.chm.lipids.pc import *
from fflip.chm.lipids.pe import *
from fflip.chm.lipids.pg import *


def parse_lipid(lipname):
    if lipname.lower() in ['dlpc', 'dppc', 'dmpc', 'popc']:
        return pc
    if lipname.lower() == 'prpc':
        return pc
    if lipname.lower() in ['dppe', 'dope', 'dlpe', 'pope']:
        return pe
    if lipname.lower()[-2:] == 'pg':
        return pg
