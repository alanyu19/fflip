#!/usr/bin/env python
#SBATCH --output=./log/%a.out
#SBATCH --error=./log/%a.err
# append additional slurm parameters here


import os
import mdtraj as md
from simtk.openmm.app import CharmmPsfFile
from rflow.trajectory import *
from fflip.analysis.scd import *


## update system info ##
# absolute paths encouraged
PSF_FILE = '../sim/plant.psf'
TRAJ_TEMPLATE = '../sim/trj/dyn{}.dcd'
FIRST_SEQNO = 101
LAST_SEQNO = 105
LAST_ATOM_INDEX_IN_UPPER_LEAFLET = 7233
SPECIAL_CARBON_WITH_SPLITTING = {
    'PLPC': ['C22'],
    'DLIPE': ['C22'],
    'PNPG': ['C22'],
    'PLPI': ['C22'],
    'DLIPS': ['C22'],
    'CET160': ['C1S', 'C5S', 'C2F'],
    'CET161': ['C1S', 'C5S', 'C2F']
}

## end of system info update ##


## start of code (no change needed below) ##
if not os.path.isdir('scd_up'):
    os.mkdir('scd_up')
if not os.path.isdir('scd_low'):
    os.mkdir('scd_low')

psf = CharmmPsfFile(PSF_FILE)
top = md.Topology.from_openmm(psf.topology)

trajs = TrajectoryIterator(
    first_sequence=FIRST_SEQNO, last_sequence=LAST_SEQNO,
    filename_template=TRAJ_TEMPLATE,
    topology_file=PSF_FILE,
    atom_selection="all", load_function=md.load_dcd
)

op_factory = OrderParameterFactory(
    top,
    special_carbons_for_splitting=SPECIAL_CARBON_WITH_SPLITTING,
    sep_leaflet=True, recipe={'hard_border': LAST_ATOM_INDEX_IN_UPPER_LEAFLET},
)

# Generate the calculation object.
# Each carbon-hydrogen vector combination is one calculation instance,
# although one calculation instance handles lipids with the same name at once.
# For example, one instance in scd_calcs may handle DLIPE C22-H2S and outputs the average.
# *** IN PRINCIPLE, USERS CAN MAKE USE OF THIS FEATURE TO ***
# *** PARALLELIZE CALCULATIONS FOR DIFFERENT LIPIDS AND CARBON POSITION ***
# *** THE EASIEST HACK WOULD BE ONLY CALLING PART OF THE scd_calcs LIST IN ONE SLURM JOB ***
scd_calcs = op_factory()

for calc in scd_calcs:
    for traj in trajs:
        calc(traj)
    hstr = ""
    for ih, h in enumerate(calc.atom2):
        if ih == 0:
            hstr += "{}".format(h)
        else:
            hstr += "_{}".format(h)
    np.savetxt(
        "./scd_up/{}-{}-{}.scd".format(calc.residue.lower(), calc.atom1.lower(), hstr.lower()), 
        np.atleast_1d(calc.scd_up)
    )
    np.savetxt(
        "./scd_low/{}-{}-{}.scd".format(calc.residue.lower(), calc.atom1.lower(), hstr.lower()), 
        np.atleast_1d(calc.scd_low)
    )
