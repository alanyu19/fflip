# -*- coding: utf-8 -*-

from fflip.chm.lipids.dppc import *
from fflip.chm.lipids.prpc import *
from fflip.chm.lipids.dmpc import *
from fflip.chm.lipids.dopc import *
from fflip.chm.lipids.pc import *
from fflip.chm.lipids.pe import *


def parse_lipid(lipname):
    if lipname.lower() in ['dlpc', 'dppc', 'dmpc', 'popc']:
        return pc
    if lipname.lower() == 'prpc':
        return pc
    if lipname.lower() in ['dppe', 'dope', 'dlpe', 'pope']:
        return pe
