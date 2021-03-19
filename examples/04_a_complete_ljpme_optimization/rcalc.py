#!/usr/bin/env python
#SBATCH --output=./.rlog/robc_%A_%a.out
#SBATCH --error=./.rlog/robc_%A_%a.err
#SBATCH --array=0-52
#SBATCH --time=01:00:00
#SBATCH --partition=sbr,ivy,hwell,k20,k40
#SBATCH --ntasks=1


import os
import sys
sys.path.append(os.getcwd())

import click
from targets import *
from fflip.ljpme.util import construct_indexes, parse_first_last


def default_first_last(target_prop):
    tp = target_prop
    if tp.prop_type == 'area' and tp.system_type == 'bilayer' and \
            tp.surface_tension == 0:
        return 51, 200
    if tp.prop_type == 'db' or tp.prop_type == 'dhh':
        return 51, 200
    if tp.prop_type == 'area' and tp.system_type == 'bilayer' and \
            tp.surface_tension != 0:
        return 61, 300
    if tp.prop_type == 'area' and tp.system_type == 'monolayer':
        return 51, 200
    if tp.prop_type == 'scd':
        return 51, 200
    if tp.prop_type == 'rdf':
        return 11, 100


@click.command()
@click.option("-i", "--property_indexes", type=str,
              help="the indexes of properties we want to calculate "
                   "(example: all / 0 / 0,2,5 / 4-10)")
@click.option("-f", "--first_trj", type=str,
              help="the first trajectory index(es)"
                   "(example: 41 / 41,61,11)")
@click.option("-l", "--last_trj", type=str,
              help="the last trajectory index(es)"
                   "(example: 100 / 100,400,200)")
@click.option("-p", "--perturbation", type=float)
@click.option("--iter", type=str, help="iteration id")
def main(property_indexes, first_trj, last_trj, perturbation, iter):

    slurm_id = os.environ["SLURM_ARRAY_TASK_ID"]

    indexes = construct_indexes(property_indexes, properties)

    if first_trj is not None and last_trj is not None:
        firsts = parse_first_last(first_trj)
        lasts = parse_first_last(last_trj)
        assert len(firsts) == len(lasts)
        assert len(firsts) <= len(indexes)
        if len(indexes) > 1 and len(firsts) == 1:
            print(
                "Broadcasting first-last trj for all targets as {}-{}".format(
                    firsts[0], lasts[0]
                )
            )
            f = firsts[0]
            l = lasts[0]
            firsts = [f for _ in range(len(indexes))]
            lasts = [l for _ in range(len(indexes))]

        elif len(indexes) > 0 and len(firsts) > 0:
            assert len(indexes) == len(firsts), \
                "Different number of starting trj and properties!"

    for count, i in enumerate(indexes):
        if int(slurm_id) == int(i):
            if first_trj is None or last_trj is None:
                first, last = default_first_last(properties[i])
            else:
                first = firsts[count]
                last = lasts[count]
            assert (last - first + 1) % 3 == 0, \
                "Please provide a trajectory range that can be divided by 3!"

            properties[i].update_percentage_of_perturbation(perturbation)
            properties[i].update_first_last_trj(first, last)
            properties[i].get_robustness(
                iteration=iter, fromfile=False,
                block_size=int((last - first + 1) / 3)
            )


if __name__ == "__main__":
    main()
