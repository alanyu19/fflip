#!/usr/bin/env python

from fflip.analysis.clustering import *
from fflip.analysis.util import find_up_low_from_crd


# UPDATE file names
# Use a (re-)centered crd!! It might be your step5_assembly.crd
PSF_FILE = '../../sim/plant.psf'
CRD_FILE = '../../sim/plant.crd'

# UPDATE the first/last trajectory indexes you want to analyze
FIRST=101
LAST=101

# UPDATE the traj template, {} replaces the index of traj
TRAJ_TEMPLATE = '../../sim/trj/dyn{}.dcd'

# UPDATE if your system contain more/different lipid types.
# The relative sizes of these and the clfix's used in ClusterLip
# can be determined by the RDFs (self and inter-species), usually
# the G-lipids are similar in size.
RADII = {
    'pe': 0.29,
    'pc': 0.29,
    'pg': 0.29,
    'pi': 0.29,
    'ps': 0.29,
    'sito': 0.29,  # sitosterol
    'cer': 0.30,  # ceramides
}   
                          

up, low = find_up_low_from_crd(CRD_FILE)

cls = ClusterLip(
    psf_file=PSF_FILE,
    # UPDATE residue names and their types
    # 'g': glycerolipids; 'c': sterols; 's': sphingo or cerimides.
    # This example is based on a sysmmetric bilayer, but you can also
    # use different top_res_info and bot_res_info for asymmetic ones.
    top_res_info=[
        ('SITO', 'c.sito'), ('PLPC', 'g.pc'), ('DLIPE', 'g.pe'),
        ('PNPG', 'g.pg'), ('PLPI', 'g.pi'), ('DLIPS', 'g.ps')
    ],
    bot_res_info=[
        ('SITO', 'c.sito'), ('PLPC', 'g.pc'), ('DLIPE', 'g.pe'),
        ('PNPG', 'g.pg'), ('PLPI', 'g.pi'), ('DLIPS', 'g.ps')
    ],
    # Don't change this recipe unless you really need to
    leaflet_recipe={'explicit_grouping': (up, low)},
    radii=RADII,
    # UPDATE min lipids to form a cluster if needed
    min_lipids=3,
    # UPDATE/DELETE/ADD pair specific if needed
    clfix={
        'sito-pc': 0.44, 'sito-pe': 0.44, 'sito-pg': 0.44, 
        'sito-pi': 0.44, 'sito-ps': 0.44,
    }
)

for i in range(FIRST, LAST + 1):
    traj = md.load_dcd(TRAJ_TEMPLATE.format(i), top=PSF_FILE)
    cls(traj, save_labels=True, first_traj_index=FIRST)

# this number can be changed to allow 
# different max in the histogram of cluster size
cls.plot_and_save(max_size=20)

cls.analyze_more(first_traj_index=FIRST, last_traj_index=LAST)
