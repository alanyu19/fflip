
from fflip.chm.lipids.dmpc import *
from fflip.chm.lipids.dppc import *
from fflip.chm.lipids.dopc import * 
from fflip.chm.lipids.prpc import * 

def parse_lipid(lipname):
    if lipname.lower() == 'dppc':
        return dppc
    if lipname.lower() == 'prpc':
        return prpc
