# -*- coding: utf-8 -*-


from fflip.drude.lipids.pc import *


def parse_lipid(lipname):
    if lipname.lower() in ['dlpc', 'dppc', 'dmpc', 'popc', 'prpc', 'pc']:
        return drude_pc
