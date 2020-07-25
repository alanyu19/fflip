# -*- coding: utf-8 -*-

#!/usr/bin/env python
#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=ivy,sbr,hwell
#SBATCH --ntasks=1


from coffe.omm.util import get_md_options as gmd

import coffe.omm.energiesforces as ef
from coffe.omm.util import *
from rflow.observables import *
from rflow.trajectory import *

option_file_name = "calc.inp"

opts = gmd(option_file_name)
starting_traj = int(opts['first'])
last_dcd = int(opts['last'])
blk_size = int(opts['blk_size']) # number of traj files each block
traj_template = str(opts['traj_template'])
psf_file = str(opts['psf'])
nlipids = int(opts['nlip'])


parameter_files = ["/u/alanyu/coffe/tests/omm/data/par_all36_lipid.prm",
                   "/u/alanyu/coffe/tests/omm/data/top_all36_lipid.rtf",
                   "/u/alanyu/coffe/tests/omm/data/toppar_water_ions.str",
                   "/u/alanyu/top_yalun/for_ua/c36ua.str" ]

psf, topology, parameters = ReadStructureParameterFiles(psf_file, parameter_files)

# The following traj is only used to set the unitcell lengths
traj = md.load_dcd(
    traj_template + "/dyn{}.dcd".format(starting_traj), top=topology
)
psf.setBox(*traj.unitcell_lengths[0])
system = psf.createSystem(
    parameters, 
    nonbondedMethod=LJPME,
    constraints=HBonds, 
    nonbondedCutoff=1.0 * u.nanometer,
    ewaldErrorTolerance=0.0005
)

trajs = TrajectoryIterator(
    first_sequence=starting_traj,
    last_sequence=starting_traj + blk_size - 1,
    filename_template = os.path.join(traj_template, 'dyn{}.dcd'),
    topology_file = psf_file,
    atom_selection = "all",
    load_function=md.load_dcd)

# this number of lipids should not be constant in the future
sa_evaluator = ef.AreaPerLipid(int(nlipids))
to_calc = TimeSeries(evaluator=sa_evaluator, filename = "./block_data/area_{}.dat".format(starting_traj))

for traj in trajs:
    to_calc(traj)

