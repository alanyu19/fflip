# -*- coding: utf-8 -*-


from fflip.drude.lipids.pc import *


def parse_lipid(lipname):
    if lipname.lower() in [
        'drude_pc', 'drude_dppc', 'drude_dmpc', 'drude_popc', 'drude_prpc'
    ]:
        return drude_pc
