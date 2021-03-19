#!/usr/bin/env python

import click
import ast
from targets import *


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

def construct_indexes(inp, properties):
    if inp=='all' or inp=='All':
        return list(range(len(properties)))
    elif ',' in inp:
        indexes = []
        values = inp.split(',')
        for v in values:
            assert 0 <= int(v) < len(properties)
            indexes.append(int(v))
        return indexes
    elif '-' in inp:
        first_last = inp.split('-')
        assert len(first_last) == 2
        first = int(first_last[0])
        last = int(first_last[1])
        assert len(properties) > last > first >= 0
        return list(range(first, last + 1))
    else:
        try:
            index = int(inp)
            assert 0 <= index < len(properties)
            return [index]
        except:
            raise Exception('Invalid property index(es)!')


def parse_first_last(inp):
    if ',' in inp:
        indexes = []
        values = inp.split(',')
        for v in values:
            indexes.append(int(v))
        return indexes
    else:
        try:
            index = int(inp)
            return [index]
        except:
            raise Exception('Invalid trajectory index(es)!')


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
@click.option(
    "-c", "--calctype", type=click.Choice(
        ['Potential', 'Observable', 'All'],case_sensitive=False
    )
)
@click.option("-p", "--perturbation", type=float,
              help='perturbation of parameter')
@click.option("--iter", type=str, help="iteration id")
@click.option("--overwrite/--keep", default=False)
@click.option("--verbose/--silent", default=True)
@click.option("-s", "--solution", default=None, help="last solution")
@click.option("-t", "--torfix", default=None, help="torsion fix file")
def main(property_indexes, first_trj, last_trj, calctype, perturbation, iter,
         overwrite, verbose, solution, torfix):
    traj_loc = '/v/gscratch/mbs/alanyu/c36ljpme'
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
        if first_trj is None or last_trj is None:
            first, last = default_first_last(properties[i])
        else:
            first = firsts[count]
            last = lasts[count]
        if verbose:
            print(
                properties[i].name,
                "(iter {}): trajectories {} to {}".format(
                    iter, first, last
                )
            )
        properties[i].update_first_last_trj(first, last)
        properties[i].update_percentage_of_perturbation(perturbation)
        if calctype is not None:
            if calctype.lower() == 'all' or calctype.lower() == 'observable':
                properties[i].calc_observable(
                    iter, traj_root=traj_loc, overwrite=overwrite
                )
            if calctype.lower() == 'all' or calctype.lower() == 'potential':
                properties[i].recalc_energy(
                    iter, traj_root=traj_loc, overwrite=overwrite, 
                    last_solution=solution, torfix=torfix
                )


if __name__ == "__main__":
    main()

