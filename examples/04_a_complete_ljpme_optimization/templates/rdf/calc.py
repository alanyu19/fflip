#!/usr/bin/env python
# -*- coding: utf-8 -*-

#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --partition=ivy,sbr,hwell
#SBATCH --ntasks=1

import sys

from rflow.trajectory import *
from coffe.omm.util import get_md_options as gmd
from fflip.analysis.rdf import RDF

option_file_name = "obscalc.inp"

opts = gmd(option_file_name)
blk_size = int(opts['block_size'])  # number of traj files each block
traj_template = os.path.join(str(opts['trj_location']), 'trj/dyn{}.dcd')
psf_file = str(opts['psf'])
lipname = str(opts['lipname']).upper()

starting_traj = int(sys.argv[1])

index = int(os.environ["SLURM_ARRAY_TASK_ID"])

psf = CharmmPsfFile(psf_file)

phos_pairs = ["Os1-HW", "O2-HW", "Os1-OW", "O2-OW"]

phos_selections = [["name O11 or name O12", "name H1 or name H2"],
                   ["name O13 or name O14", "name H1 or name H2"],
                   ["name O11 or name O12", "name OH2"],
                   ["name O13 or name O14", "name OH2"]]

es_pairs = ["Ob-OW", "Os2-OW", "Ob-HW", "Os2-HW"]

es_selections = [["name O22 or name O32", "name OH2"],
                 ["name O21 or name O31", "name OH2"],
                 ["name O22 or name O32", "name H1 or name H2"],
                 ["name O21 or name O31", "name H1 or name H2"]]

pairs = phos_pairs + es_pairs
selections = phos_selections + es_selections


def get_rdf_mdtraj(
        rdf_name, rdf_atom_selection, rdf_psf, trajectories,
        start_traj=starting_traj
):
    rdf_mdtraj = RDF(
        rdf_psf, rdf_name, rdf_atom_selection, r_range=[0, 1.0], dimension=3,
        method='mdtraj'
    )
    for trj in trajectories:
        rdf_mdtraj(trj, start_traj, save_sparse=True, sparse_bin_width=0.003,
                   sparse_subfolder='./block_data')
    rdf_mdtraj.save_to_file()


trajs = TrajectoryIterator(
    first_sequence=starting_traj, last_sequence=starting_traj,
    # this is currently taking only one trj
    filename_template=traj_template,
    topology_file=psf_file,
    atom_selection="all",
    load_function=md.load_dcd
)

name = pairs[index-1]
atom_selection = selections[index-1]
get_rdf_mdtraj(name, atom_selection, psf, trajs)
