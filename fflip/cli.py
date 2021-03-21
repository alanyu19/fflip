# -*- coding: utf-8 -*-

"""Console script for fflow."""

import os
import sys
sys.path.append(os.getcwd())

import click
from targets import *
import fflip
from fflip.ljpme.util import construct_indexes, parse_first_last


@click.group()
@click.version_option(version=fflip.__version__)
@click.pass_context
def main(ctx, args=None):
    """Console script for fflip"""
    click.echo("FFLiP: a python package for CHARMM/Drude lipid force field parameterization")
    return 0


@main.command()
@click.option("-i", "--property_indexes", type=str,
              help="the indexes of properties we want to calculate "
                   "(example: all / 0 / 0,2,5 / 4-10)")
@click.option("-f", "--first_trj", type=str,
              help="the first trajectory index(es)"
                   "(example: 41 / 41,61,11)")
@click.option("-l", "--last_trj", type=str,
              help="the last trajectory index(es)"
                   "(example: 100 / 100,400,200)")
@click.option("-p", "--perturbation", type=float,
              help='the perturbation size used in reweighting (%)')
@click.option("--iteration", type=str, help="iteration id")
@click.option("-P", "--partition", type=str, help="slurmq partition")
def scalc(property_indexes, first_trj, last_trj, perturbation, iteration, partition):
    """
    Perform the sensitivity calculation
    """
    
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

    from threading import Thread
    for count, i in enumerate(indexes):
        if first_trj is None or last_trj is None:
            first, last = properties[i].option_scheme.critical_seqno
        else:
            first = firsts[count]
            last = lasts[count]
        properties[i].update_percentage_of_perturbation(perturbation)
        properties[i].update_first_last_trj(first, last)
        thread = Thread(
            target=properties[i].get_sensitivity,
            kwargs={'iteration': iteration, 'force_redo': True, 'partition': partition}
        )
        thread.start()
        # properties[i].get_sensitivity(iteration, force_redo=True)


@main.command()
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
@click.option("--iteration", type=str, help="iteration id")
@click.option("-P", "--partition", type=str, help="slurmq partition")
def rcalc(property_indexes, first_trj, last_trj, perturbation, iteration, partition):
    """
    Perform the error(robustness) analysis on sensitivities
    """
    
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

    from threading import Thread
    for count, i in enumerate(indexes):
        if first_trj is None or last_trj is None:
            first, last = default_first_last(properties[i])
        else:
            first = firsts[count]
            last = lasts[count]
        assert (last - first + 1) % 3 == 0, \
            "Please provide a trajectory range that can be divided by 3!"

        properties[i].update_percentage_of_perturbation(perturbation)
        properties[i].update_first_last_trj(first, last)
        thread = Thread(
            target=properties[i].get_robustness,
            kwargs={'iteration': iteration, 
                    'fromfile': False,
                    'block_size': int((last - first + 1) / 3),
                    'partition': partition}
        )
        thread.start()


@main.command()
@click.option("-l", "--location", type=str,
              help="root directory to run the simulation")
@click.option("-s", "--sfile", type=str, default=None,
              help="solution file for nonbonded parameters in 1d array format")
@click.option("-t", "--tfile", type=str, default=None,
              help="file for torsion fixes")
@click.option("-i", "--iteration", type=str,
              help="interation id, can be 0, 1, 2 ... or any string")
def simulate(location, sfile, tfile, iteration):
    """
    Initiate the simulations 
    """
    # from targets import *
    for prop in properties:
        prop.simulate(
            iteration=iteration, trj_folder=location,
            # change change_para is set to True,
            # make sure to provide a solution file (sfile)
            change_para=False, torfix_file=tfile, solution_file=sfile
        )


@main.command()
@click.option("-i", "--property_indexes", type=str,
              help="the indexes of properties we want to calculate "
                   "(example: all / 0 / 0,2,5 / 4-10)")
@click.option("-tl", "--traj_loc", type=str,
              help="location of the trajectory/simulation")
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
@click.option("--iteration", type=str, help="iteration id")
@click.option("--overwrite/--keep", default=False)
@click.option("--verbose/--silent", default=True)
@click.option("-s", "--solution", default=None, help="last solution")
@click.option("-t", "--torfix", default=None, help="torsion fix file")
def obspot(property_indexes, traj_loc, first_trj, last_trj, calctype, perturbation, iteration,
           overwrite, verbose, solution, torfix):
    assert traj_loc is not None
    traj_loc = os.path.abspath(traj_loc)
    assert property_indexes is not None
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
            first, last = properties[i].option_scheme.critical_seqno
        else:
            first = firsts[count]
            last = lasts[count]
        if verbose:
            print(
                properties[i].name,
                "(iteration {}): trajectories {} to {}".format(
                    iteration, first, last
                )
            )
        properties[i].update_first_last_trj(first, last)
        properties[i].update_percentage_of_perturbation(perturbation)
        if calctype is not None:
            if calctype.lower() == 'all' or calctype.lower() == 'observable':
                properties[i].calc_observable(
                    iteration, traj_root=traj_loc, overwrite=overwrite
                )
            if calctype.lower() == 'all' or calctype.lower() == 'potential':
                properties[i].recalc_energy(
                    iteration, traj_root=traj_loc, overwrite=overwrite, 
                    last_solution=solution, torfix=torfix
                )

@main.command()
@click.option("-i", "--iteration", type=str,
              help="iteration index/label")
@click.option("-p", "--perturbation", type=str,
              help="perturbation used in the reweighting")
@click.option("-s", "--sigrst", type=float, default=0.05,
              help="restraint on sigma, default is 0.05, increase for more restraint")
@click.option("-e", "--epsrst", type=float, default=0.05,
              help="restraint on epsilon, defult is 0.05, increase for more restraint")
@click.option("-c", "--chrgrst", type=float, default=0.025,
              help="restraint on charge, defualt is 0.025, increase for more restraint")
@click.option("-u", "--uncertainty_scaling", type=str, default=500,
              help="increase to apply more restraint")
@click.option("--ssr", is_flag=True, help="plot properties' contributions to SSR")
def linearopt(iteration, perturbation, sigrst, epsrst, chrgrst, uncertainty_scaling, ssr):
    print('there are {} normal properties and {} special properties'.format(
        len(properties), len(special_properties))
    )

    le = PropertyLinearEstimator(
        properties, 
        special_properties, 
        uncertainty_scaling=uncertainty_scaling  # larger scaling should be used for more restraint
    )
    
    for prop in le.target_properties:
        prop.update_percentage_of_perturbation(perturbation)
        prop.get_sensitivity(iteration, use_cluster=False, quit=True)
    for prop in le.special_properties:
        prop.get_sensitivity(iteration)
        prop.update_scaling_using_exp()
    
    for prop in le.target_properties:
        prop.get_robustness(iteration)
    
    # This should be an option (TODO: Yalun)
    le.get_qm_charge_from_file(root + '/qm/qm.csv')
    le.get_sensitivity_matrix()
    le.get_weight_matrix(
        qm_weight=0.00, hard_bounds={
            'sigma': sigrst, 'epsilon': epsrst, 'charge': chrgrst
        },
        # This is currently hard-coded
        drop_bounds={'sigma': 0.2, 'epsilon': 0.2, 'charge': 0.2},
    ) 
    le.get_deviation_vector()

    # np.savetxt('W.mtx', le.W)
    # np.savetxt('S.mtx', le.S)
    # np.savetxt('F.mtx', le.F)
    
    solution = le(save_result=False)
    le.update_weight(
        hard_bounds={'sigma': sigrst, 'epsilon': epsrst, 'charge': chrgrst},
        # parameter change smaller than this would be dropped
        lower_bound=0.2
    )
    if not os.path.isdir("solutions"):  # directory for slurm output files
        os.mkdir("solutions")
    # TODO: this should be enhanced!
    solution = le(save_result=ssr, ssr_file = './solutions/ssr_' + iteration + '.png')
    np.savetxt('./solutions/solution_{}.txt'.format(iteration), solution)

    # TODO (yalun): this should be supported!
    # param = dppc.parse_nbgroups(groups=properties[0].parse_groups)
    # previous = np.loadtxt('./solutions/c2aa/solution.txt')
    # final = combine_solutions([previous, solution], param)
    # np.savetxt('solution.txt', final)


def entrypoint():
    sys.exit(main(obj={}))

