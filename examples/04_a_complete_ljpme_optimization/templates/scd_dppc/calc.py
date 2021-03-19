#!/usr/bin/env python

# -*- coding: utf-8 -*-

#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=ivy,sbr,hwell
#SBATCH --ntasks=1

from simtk.openmm.app import LJPME

from coffe.omm.util import get_md_options as gmd
from coffe.omm.util import *
from rflow.observables import *
from rflow.trajectory import *
from fflip.analysis.scd import OrderParameterCalculation as opc


option_file_name = "obscalc.inp"

opts = gmd(option_file_name)
blk_size = int(opts['block_size'])  # number of traj files each block
traj_template = os.path.join(str(opts['trj_location']), 'trj/dyn{}.dcd')
psf_file = str(opts['psf'])
lipname = str(opts['lipname']).upper()

starting_traj = int(os.environ['SLURM_ARRAY_TASK_ID'])


parameter_files = ["/u/alanyu/coffe/tests/omm/data/par_all36_lipid.prm",
                   "/u/alanyu/coffe/tests/omm/data/top_all36_lipid.rtf",
                   "/u/alanyu/coffe/tests/omm/data/toppar_water_ions.str"]

psf, topology, parameters = read_structure_parameter_files(
    psf_file, parameter_files
)

traj = md.load_dcd(
    traj_template.format(starting_traj), top=topology
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
    filename_template=traj_template,
    topology_file=psf_file,
    atom_selection="all",
    load_function=md.load_dcd
)

scd_evaluators = []
for i in range(3, 16):
    scd_evaluators.append(
        opc(
            topology, lipname, 'C2{}'.format(i),
            ['H{}R'.format(i), 'H{}S'.format(i)]
        )
    )
    scd_evaluators.append(
        opc(
            topology, lipname, 'C3{}'.format(i),
            ['H{}X'.format(i), 'H{}Y'.format(i)]
        )
    )

scd_evaluators_extra = list()
scd_evaluators_extra.append(opc(topology, lipname, 'C22', ['H2R']))
scd_evaluators_extra.append(opc(topology, lipname, 'C22', ['H2S']))
scd_evaluators.append(opc(topology, lipname, 'C32', ['H2X', 'H2Y']))

# Headgroup (compare to Toward Atomistic Resolution Structure of
# Phosphoatidylcholine Headgroup and Glycerol Backbone at Different
# Ambient Conditions. JPCB 2015, 119, 15075-15088)
scd_evaluators_head = list()
scd_evaluators_head.append(opc(topology, lipname, 'C11', ['H11A', 'H11B']))
scd_evaluators_head.append(opc(topology, lipname, 'C12', ['H12A', 'H12B']))
# Glycerol Backbone
scd_evaluators_gly = list()
scd_evaluators_gly.append(opc(topology, lipname, 'C1', ['HA']))
scd_evaluators_gly.append(opc(topology, lipname, 'C1', ['HB']))
scd_evaluators_gly.append(opc(topology, lipname, 'C2', ['HS']))
scd_evaluators_gly.append(opc(topology, lipname, 'C3', ['HX']))
scd_evaluators_gly.append(opc(topology, lipname, 'C3', ['HY']))

scd = []
for scd_evaluator in scd_evaluators:
    scd.append(
        TimeSeries(
            evaluator=scd_evaluator,
            filename="./block_data/scd-{}_{}.dat".format(
                scd_evaluator.atom1.lower(), starting_traj
            )
        )
    )

for scd_evaluator in scd_evaluators_extra:
    scd.append(
        TimeSeries(
            evaluator=scd_evaluator,
            filename="./block_data/scd-{}-{}_{}.dat".format(
                scd_evaluator.atom1.lower(), scd_evaluator.atom2[0].lower(),
                starting_traj
            )
        )
    )

for scd_evaluator in scd_evaluators_head:
    scd.append(
        TimeSeries(
            evaluator=scd_evaluator,
            filename="./block_data/scd-{}_{}.dat".format(
                scd_evaluator.atom1.lower(), starting_traj
            )
        )
    )

for scd_evaluator in scd_evaluators_gly:
    scd.append(
        TimeSeries(
            evaluator=scd_evaluator,
            filename="./block_data/scd-{}-{}_{}.dat".format(
                scd_evaluator.atom1.lower(), scd_evaluator.atom2[0].lower(),
                starting_traj
            )
        )
    )

for traj in trajs:
    for s in scd:
        s(traj)

