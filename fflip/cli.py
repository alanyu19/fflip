# -*- coding: utf-8 -*-

"""Console script for fflow."""

import os
import sys
sys.path.append(os.getcwd())

import click
if os.path.isfile('targets.py'):
    from targets import *
if os.path.isfile('models.py'):
    from models import *
import fflip
from fflip.ljpme.util import construct_indexes, parse_first_last, dict_to_str
# analysis related:
import mdtraj as md
import rflow.observables as obs
from rflow.trajectory import TrajectoryIterator
from fflip.plot.area_plot import plot_area
import fflip.analysis.area_analysis as aa


@click.group()
@click.version_option(version=fflip.__version__)
@click.pass_context
def main(ctx, args=None):
    """Console script for fflip"""
    click.echo(
        "FFLiP: A Python Package for CHARMM/Drude Lipid FF Parameterization"
    )
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
@click.option("-bsp", "--block_size_potential", type=int, default=-1,
              help="block size (# traj) for potential calculation")
@click.option("-bso", "--block_size_observable", type=int, default=-1,
              help="block size (# traj) for observable calculation")
@click.option("-p", "--perturbation", type=float,
              help='the perturbation size used in reweighting (%)')
@click.option("--iteration", type=str, help="iteration id")
@click.option("-P", "--partition", type=str, help="slurmq partition")
def scalc(property_indexes, first_trj, last_trj, block_size_potential,
          block_size_observable, perturbation, iteration, partition):
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
        properties[i].update_amount_of_perturbation(perturbation)
        properties[i].update_first_last_trj(first, last)
        if block_size_potential > 0:
            properties[i].update_block_sizes(potential=block_size_potential)
        if block_size_observable > 0:
            properties[i].update_block_sizes(observable=block_size_observable)
        if count % 3 == 0:
            if count != 0:
                for thread_ in threads:
                    thread_.join()
            threads = []
        thread = Thread(
            target=properties[i].get_sensitivity,
            kwargs={'iteration': iteration, 'force_redo': True, 'partition': partition}
        )
        thread.start()
        threads.append(thread)
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
@click.option("-bsp", "--block_size_potential", type=int, default=-1,
              help="block size (# traj) for potential calculation")
@click.option("-bso", "--block_size_observable", type=int, default=-1,
              help="block size (# traj) for observable calculation")
@click.option("-p", "--perturbation", type=float)
@click.option("--iteration", type=str, help="iteration id")
@click.option("-P", "--partition", type=str, help="slurmq partition")
def rcalc(property_indexes, first_trj, last_trj, block_size_potential,
          block_size_observable, perturbation, iteration, partition):
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

        properties[i].update_amount_of_perturbation(perturbation)
        properties[i].update_first_last_trj(first, last)
        if block_size_potential > 0:
            properties[i].update_block_sizes(potential=block_size_potential)
        if block_size_observable > 0:
            properties[i].update_block_sizes(observable=block_size_observable)
        if count % 3 == 0:
            if count != 0:
                for thread_ in threads:
                    thread_.join()
            threads = []
        thread = Thread(
            target=properties[i].get_robustness,
            kwargs={'iteration': iteration, 
                    'fromfile': False,
                    'block_size': int((last - first + 1) / 3),
                    'partition': partition}
        )
        thread.start()
        threads.append(thread)
        # properties[i].get_robustness(
        #     iteration=iteration, fromfile=False, partition=partition,
        #     block_size = int((last - first + 1) / 3)
        # )


@main.command()
@click.option("-l", "--location", type=str,
              help="root directory to run the simulations")
@click.option("-t", "--time", type=int,
              help="simulaiton length (in sequence number, real time depend on run.py)")
@click.option("-c", "--change_parameters", default='no', type=click.Choice(['yes', 'no']),
              help="whether to change nonbonded parameters, if yes, sfile is needed")
@click.option("-sf", "--sfile", type=str, default=None,
              help="solution file for nonbonded parameters in 1d array format")
@click.option("-tf", "--tfile", type=str, default=None,
              help="file for torsion fixes")
@click.option("-x", "--boxx", type=float, help="simulation box size in x/y (Å)")
@click.option("-z", "--boxz", type=float, help="simulation box size in z (Å)")
@click.option("--zmode", type=int, default=0,
              help="Z mode for OpenMM MC Membrane Barostat")
@click.option("-i", "integrator", type=click.Choice(['L', 'N']), default='L',
              help="Integrator: Langevin (L) or Nose-Hoover (N)")
@click.option("-b", "--barostat", type=click.Choice(['MCM', 'MC']), default='MCM',
              help="Integrator: MC Membrane (MCM) or MonteCarlo (MC)")
@click.option("-I", "--iteration", type=str,
              help="interation id, can be 0, 1, 2 ... or any string")
@click.option("--start", is_flag=True,
              help="start immediately? (otherwise go to sim folders to start manually)")
@click.option("-v", "--verbose", type=int, default=1)
@click.option("--overwrite", is_flag=True)
def simulate(location, time, change_parameters, sfile, tfile,
             boxx, boxz, zmode, integrator, barostat, iteration,
             start, verbose, overwrite):
    """
    Initiate the simulations 
    """
    cp_match = {'yes': True, 'no': False}
    for prop in properties:
        prop.simulate(
            iteration=iteration, trj_folder=location, last_seqno=time,
            change_param=cp_match[change_parameters],
            torfix_file=tfile, solution_file=sfile,
            boxx=boxx, boxz=boxz, zmode=zmode,
            barostat=barostat, integrator=integrator, start=start,
            verbose=verbose, overwrite=overwrite
        )


@main.command()
@click.option("-l", "--location", type=str,
              help="root directory to run the simulations")
@click.option("-c", "--change_parameters", default='no', type=click.Choice(['yes', 'no']),
              help="whether to change nonbonded parameters, if yes, sfile is needed")
@click.option("-sf", "--sfile", type=str, default=None,
              help="solution file for nonbonded parameters in 1d array format")
@click.option("-tf", "--tfile", type=str, default=None,
              help="file for torsion fixes")
@click.option("-I", "--iteration", type=str,
              help="interation id, can be 0, 1, 2 ... or any string")
@click.option("--toppar", type=str, default=None, help="toppar path")
@click.option("--start", is_flag=True,
              help="start immediately? (otherwise go to sim folders to start manually)")
def simulate_mc(location, change_parameters, sfile, tfile, iteration, toppar, start):
    """
    Initiate the simulations 
    """
    cp_match = {'yes': True, 'no': False}
    for mcp in mcps:
        mcp.generate_model_compound()
        trj_folder = os.path.join(location, 'iter{}'.format(iteration), 'model_compound', mcp.name.lower())
        # TODO: add other parameters
        mcp.simulate(trj_folder, toppar, start)


@main.command()
@click.option("-l", "--location", type=str,
              help="root directory of the simulations")
@click.option("-I", "--iteration", type=str,
              help="interation id, can be 0, 1, 2 ... or any string")
def dihedral_mc(location, iteration):
    """
    Dihedral Distribution
    """
    for mcp in mcps:
        mcp.generate_model_compound()
        trj_folder = os.path.join(location, 'iter{}'.format(iteration), 'model_compound', mcp.name.lower())
        mcp.dihedral_distribution(trj_folder)


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
@click.option("-bsp", "--block_size_potential", type=int, default=-1,
              help="block size (# traj) for potential calculation")
@click.option("-bso", "--block_size_observable", type=int, default=-1,
              help="block size (# traj) for observable calculation")
@click.option(
    "-c", "--calctype", type=click.Choice(
        ['Potential', 'Observable', 'All'],case_sensitive=False
    )
)
@click.option("-p", "--perturbation", type=float,
              help='perturbation of parameter (suggestion: 0.001)')
@click.option("--iteration", type=str, help="iteration id")
@click.option("--overwrite/--keep", default=False)
@click.option("--verbose/--silent", default=True)
@click.option("--toppar", type=str, default=None, help="toppar path")
@click.option("-s", "--solution", default=None, help="last solution")
@click.option("-t", "--torfix", default=None, help="torsion fix file")
def obspot(property_indexes, traj_loc, first_trj, last_trj, block_size_potential,
           block_size_observable, calctype, perturbation, iteration, overwrite,
           verbose, toppar, solution, torfix):
    """
    Calculate energies (original+perturbed) and observables
    """
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
        properties[i].update_amount_of_perturbation(perturbation)
        if block_size_potential > 0:
            properties[i].update_block_sizes(potential=block_size_potential)
        if block_size_observable > 0:
            properties[i].update_block_sizes(observable=block_size_observable)
        if calctype is not None:
            if calctype.lower() == 'all' or calctype.lower() == 'observable':
                properties[i].calc_observable(
                    iteration, traj_root=traj_loc, overwrite=overwrite,
                )
            if calctype.lower() == 'all' or calctype.lower() == 'potential':
                toppar = os.path.realpath(toppar)
                properties[i].recalc_energy(
                    iteration, traj_root=traj_loc, overwrite=overwrite, 
                    last_solution=solution, torfix=torfix,
                    toppar_path=toppar
                )

@main.command()
@click.option("-i", "--iteration", type=str,
              help="iteration index/label")
@click.option("-p", "--perturbation", type=str,
              help="perturbation used in the reweighting")
@click.option("-s", "--sigrst", type=float, default=5,
              help="restraint on sigma, default is 5, increase for more restraint")
@click.option("-e", "--epsrst", type=float, default=5,
              help="restraint on epsilon, defult is 5, increase for more restraint")
@click.option("-t", "--tlrst", type=float, default=4,
              help="restraint on thole, default is 4, increase for more restraint")
@click.option("-a", "--aprst", type=float, default=4,
              help="restraint on alpha, defult is 4, increase for more restraint")
@click.option("-c", "--chrgrst", type=float, default=3,
              help="restraint on charge, defualt is 3, increase for more restraint")
@click.option("-u", "--uncertainty_scaling", type=float, default=500,
              help="default is 500, increase to apply more restraint")
@click.option("--hasqm", is_flag=True,
              help="if QM partial charges are available ($root/qm/qm.csv)")
@click.option("--qmw", type=float, default=None,
              help="weight of qm charges, performance not tested!")
@click.option("--previous", type=str, default=None, help="previous iteration index/label")
@click.option("--ssr", is_flag=True, help="plot properties' contributions to SSR")
def linearopt(iteration, perturbation, sigrst, epsrst, chrgrst, tlrst, aprst, uncertainty_scaling, hasqm, qmw, previous, ssr):
    print('there are {} normal properties and {} special properties'.format(
        len(properties), len(special_properties))
    )

    le = PropertyLinearEstimator(
        properties, 
        special_properties, 
        uncertainty_scaling=uncertainty_scaling  # larger scaling should be used for more restraint
    )
    # le.exchange_exp_value('scd-c22-h2r', 'scd-c22-h2s')
    
    for prop in le.target_properties:
        prop.update_amount_of_perturbation(perturbation)
        prop.get_sensitivity(iteration, use_cluster=False, quit=True)
    for prop in le.special_properties:
        prop.get_sensitivity(iteration)
        prop.update_scaling_using_exp()
    
    for prop in le.target_properties:
        prop.get_robustness(iteration)
    
    # This should be an option (TODO: Yalun)
    if hasqm:
        le.get_qm_charge_from_file(root + '/qm/qm.csv')
        assert qmw is not None
    # else:
    #     le.zero_qm_charge()
    #     qmw = 0
    le.gen_sensitivity_matrix()
    le.gen_weight_matrix(
        hard_bounds={
            'sigma': sigrst, 'epsilon': epsrst, 'charge': chrgrst,
            'thole': tlrst, 'alpha': aprst
        },
        # This is currently hard-coded
        drop_bounds={
            'sigma': 1, 'epsilon': 1, 'charge': 1,
            'alpha': 1, 'thole': 1},
        qmc_weight=qmw, verbose=True,
        # forbid={'charge': ['C2', 'C3', 'H2R', 'H2S', 'H2X', 'H2Y']}
    ) 
    le.gen_target_vector()
    # mkdir for result
    if not os.path.isdir("solutions"):  # directory for slurm output files
        os.mkdir("solutions")
    if not os.path.isdir("solutions/iter{}".format(iteration)):  # directory for slurm output files
        os.mkdir("solutions/iter{}".format(iteration))
    result_path = "solutions/iter{}".format(iteration)
    # save matrixes
    np.savetxt(os.path.join(result_path, "w.mtx"), le.W)
    np.savetxt(os.path.join(result_path, "s.mtx"), le.S)
    np.savetxt(os.path.join(result_path, "t.mtx"), le.T)
    solution = le(save_result=False)
    le.update_weight(
        hard_bounds={'sigma': sigrst, 'epsilon': epsrst, 'charge': chrgrst,
                     'thole': tlrst, 'alpha': aprst},
        # parameter change (measured in % or 0.01 e) 
        # smaller than this would be dropped
        min_change=0.002
    )
    solution = le(
        save_result=ssr, ssr_file=os.path.join(result_path, "ssr_" + iteration + ".png"),
        result_file=os.path.join(result_path, "predicted_{}.csv".format(iteration))
    )
    # save essential optimization parameters
    with open(
        os.path.join(result_path, "optimization_info.txt"), "w" 
    ) as fw:
        fw.write("uncertainty_scaling:\n")
        fw.write(str(le.uncertainty_scaling) + "\n")
        fw.write("----------------\n")
        fw.write("hard_bounds:\n")
        fw.write(dict_to_str(le.hard_bounds))
        fw.write("----------------\n")
        fw.write("drop_bounds:\n")
        fw.write(dict_to_str(le.drop_bounds))
        fw.write("----------------\n")
        fw.write("forbidden:\n")
        fw.write(dict_to_str(le.forbid))
        fw.write("----------------\n")
        fw.write("min_change:\n")
        fw.write(str(le.min_change))
    # write out solutions
    if previous is None:
        np.savetxt(os.path.join(result_path, "solution_{}.txt".format(iteration)), solution)
        with open(
            os.path.join(result_path, "solution_{}_with_param.txt".format(iteration)), 'w'
        ) as fw:
            for param, s in zip(le.parameters, solution):
                fw.write("{0:5s} {1:10s} {2:>9.4f}\n".format(param.center_names[0], param.par_type, s))
    else:
        # currently only support one lipid type at once
        lipid_ = le.target_properties[0].lipid
        nbgroups_ = le.target_properties[0].groups
        param = lipid_.parse_nbgroups(groups=nbgroups_)
        previous = np.loadtxt('./solutions/solution_{}.txt'.format(previous))
        final = combine_solutions([previous, solution], param)
        np.savetxt('./solutions/solution_{}.txt'.format(iteration), final)
        with open(
            os.path.join(result_path, "solution_{}_with_param.txt".format(iteration)), 'w'
        ) as fw:
            for param, s in zip(le.parameters, final):
                fw.write("{0:5s} {1:10s} {2:>9.4f}\n".format(param.center_names[0], param.par_type, s))


@main.command()
@click.option('-p', '--psf_file', type=str, help='PSF file used to run the simulation')
@click.option('-l', '--nlip', default=36, help='Number of lipids. Default is 36.')
@click.option('-t', '--traj_template', default='trj/dyn{}.dcd',
    help='Trajetory file template, use curly brace for index, default is "trj/dyn{}.dcd"')
# @click.option('-b', '--block_size', default=10, type=float)
@click.option('-f', '--first_seqno', default=1, type=int, help='first seqno to calculate the area')
@click.option('-e', '--last_seqno', default=None, type=int, help='last seqno to calculate the area (optional)')
def area_per_lipid(psf_file, nlip, traj_template, first_seqno, last_seqno=None):
    """ Area per Lipid """
    root_dir = os.path.dirname(os.path.realpath(psf_file))
    print(root_dir)
    if last_seqno is None:
        next_seqno_file = os.path.join(root_dir, "next.seqno")
        print(next_seqno_file)
        assert os.path.isfile(next_seqno_file), "No `next.seqno` file detected, provide last_seqno instead!"
        with open(next_seqno_file) as f:
            next_seqno = int(f.readline())
            last_seqno = int(next_seqno - 1)
    trajs = TrajectoryIterator(
        first_sequence=first_seqno,
        last_sequence=last_seqno,
        filename_template=traj_template,
        topology_file=psf_file,
        atom_selection="all",
        load_function=md.load_dcd
    )

    sa_evaluator = obs.AreaPerLipid(int(nlip))
    to_calc = obs.TimeSeries(
        evaluator=sa_evaluator, filename="area.dat"
    )

    for traj in trajs:
        to_calc(traj)


@main.command()
@click.option("-i", "--input_file", type=str, help='Input area data file')
@click.option("-l", "--nlip", type=int, help='Number of lipids per leaflet')
@click.option("-t", "--temperature", type=float, help='Temperature of the simulation (Kelvin)')
@click.option("-b", "--block_size", type=float, default=10, help='Block size to use for statistical analysis,\
    for phospholipids 10 ns is recommended (also is default)')
@click.option("-s", "--interval", type=float, help='Time interval (ns) between data points')
@click.option("-f", "--force", is_flag=True, help='If this flag provided, force the calculation even if not sufficient equlibrium data')
@click.option("-ff", "--force_fraction", type=float, default=0.6, help='Proportion of data forced to use (ignoring the equilibrium starting point)')
def analyze_area_per_lipid(input_file, nlip, temperature, block_size, interval, force, force_fraction):
    """
    Equilibrium detection and area per lipid plot using area per lipid time series data
    """
    aa.get_area_and_ka(nlipid=nlip, temperature=temperature, block_size=block_size, step_size=interval,
                       force_calc=force, force_fraction=force_fraction, file_to_read=input_file)


def entrypoint():
    sys.exit(main(obj={}))
