# -*- coding: utf-8 -*-

from fflip.chm.lipids.dppc import *
from fflip.chm.lipids.prpc import *
from fflip.chm.lipids.dmpc import *
from fflip.chm.lipids.dopc import *
from fflip.chm.lipids.pe import *


def parse_lipid(lipname):
    if lipname.lower() in ['dlpc', 'dppc', 'dmpc']:
        return dppc
    if lipname.lower() == 'prpc':
        return prpc
    if lipname.lower() in ['dppe', 'dope', 'dlpe', 'pope']:
        return pe
