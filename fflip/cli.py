# -*- coding: utf-8 -*-

"""Console script for rickflow."""

import os
import sys
import click

#import rflow
#from rflow import TrajectoryIterator, Distribution, TransitionCounter, make_topology, center_of_mass_of_selection
#import mdtraj as md
import numpy as np

@click.group()
@click.version_option(version=fflip.__version__)
def main(args=None):
    """Console script for rickflow."""
    click.echo("Rickflow: a python package to facilitate running and analyzing jobs in OpenMM using CHARMM defaults.")
    return 0


@main.command()
@click.argument("batch", type=click.Path(exists=True))
def submit(batch):
    """
    Submit the workflow using a batch script.
    """
    if not os.path.isdir("log"):  # directory for slurm output files
        os.mkdir("log")
    assert os.path.isfile(batch)
    cwd = os.path.basename(os.getcwd())
    os.system("sbatch -o log/{}-%j.log {}".format(cwd, batch))



if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
