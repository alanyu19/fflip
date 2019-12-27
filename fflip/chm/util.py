
from fflip.chm.lipids.dmpc import *
from fflip.chm.lipids.dppc import *
from fflip.chm.lipids.dopc import * 
from fflip.chm.lipids.prpc import * 

def parse_lipid(lipname):
    if lipname.lower() in ['dlpc', 'dppc', 'dmpc']:
        return dppc
    if lipname.lower() == 'prpc':
        return prpc
