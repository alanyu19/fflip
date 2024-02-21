#!/usr/bin/env python

from fflip.analysis.clustering import *
from fflip.analysis.util import find_up_low_from_crd


# UPDATE PSF/CRD, the *.crd file should be (re)centered
PSF_FILE = '../../sim/plant.psf'
CRD_FILE = '../../sim/plant.crd'

# UPDATE the first/last trajectory indexes you want to analyze
FIRST=101
LAST=101
BLOCKSIZE=1

# UPDATE the traj template, {} replaces the index of traj
TRAJ_TEMPLATE = '../../sim/trj/dyn{}.dcd'
NFRAME_PER_TRAJ = 100

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

# Change the inverval (5) to have sparser/denser output.
# 5 here means save one snapshot (for both upper and lower)
# for every 5 trajectories (the last traj gets used).
# Note the "frame" can be set to any number smaller than
# the frame number in each trajectory, though it is usually
# set to the last frame.
# "edge" controls the range of the plot (please use 1~3).
# 1.5 means the edge of the plot area is 1.5 * sim.box.edge.
# "addlabel" controls x/y labeling, you might want to remove labels
# (set to False) if you want to combine plots and add shared labels
# later on in some graphic software.
# "highcontrast" controls the color of lipids not in a frame, use "True"
# leads to grey used for these lipids. "False" leads to same coloring as
# those in cluster but non-fill circles.
# * Colors of different lipid types are controlled 
# in clustering.py (first function in that code).
for i in range(FIRST, LAST + 1, BLOCKSIZE):
    traj = md.load_dcd(TRAJ_TEMPLATE.format(i), top=PSF_FILE)
    cls.write_xy(traj, itraj=i+BLOCKSIZE-1, frame=NFRAME_PER_TRAJ)
    cls.plot_cluster(itraj=i+BLOCKSIZE-1, frame=NFRAME_PER_TRAJ, edge=1.5, addlabel=True, highcontrast=False)

