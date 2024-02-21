import os                                   
import numpy as np
from simtk.openmm.app import CharmmPsfFile 
import mdtraj as md
from fflip.analysis.rdf import *


PSF_FILE = "../../plant.psf"
CRD_FILE = "../../plant.crd"
FIRST = 101
LAST = 101
TRAJ_TEMPLATE = '../../sim/trj/dyn{}.dcd'
# LAST_ATOM_INDEX_UPPER_LEAFLET = 99999 # needed when using Method 1 below.
# change blocksize (number of trajs per block) if needed
BLOCKSIZE = 1

psf = CharmmPsfFile(PSF_FILE)

sito_plpc_rdf = RDF(
    psf, 'plpc_sito', ["resname SITO and name O3", "resname PLPC and name P"], 
    r_range=(0, 2), dimension=2, separate_leaflet=True
)

# Method 1. using user-given last atom index of upper leaflet
# for i in range(FIRST, LAST + 1):
#     trj = md.load_dcd(TRAJ_TEMPLATE.format(i), top=PSF_FILE)
#     sito_plpc_rdf(
#         trj,
#         first_trj_index=FIRST,
#         up_low_recipe={'hard_border': LAST_ATOM_INDEX_UPPER_LEAFLET}, 
#         # the number above is found from psf for separating upper & lower,
#         # doesn't work well with membrane built from CHARMM-GUI containing
#         # Glycolipids (maybe more restrictions).
#         save_blocks_interval=BLOCKSIZE,
#         save_to='./data'
#     )


# Method 2. detecting upper/lower leaflet atoms automatically
up, low = find_up_low_from_crd(CRD_FILE)
for i in range(FIRST, LAST + 1):
    trj = md.load_dcd(TRAJ_TEMPLATE.format(i), top=PSF_FILE)
    sito_plpc_rdf(
        trj,
        first_trj_index=FIRST,
        up_low_recipe={'explicit_grouping': (up, low)}, 
        # the number above is found from psf for separating upper & lower,
        # doesn't work well with membrane built from CHARMM-GUI containing
        # Glycolipids (maybe more restrictions).
        save_blocks_interval=BLOCKSIZE,
        save_to='./data'
    )
