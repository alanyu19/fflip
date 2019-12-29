# -*- coding: utf-8 -*-

import scipy.optimize as sopt
from coffe.omm.torsionfuncs import *


class ObjfuncDihedral(object):
    def __init__(self, sim, ref, reweighter):
        self.sim_data = sim
        self.ref_data = ref
        self.reweighter = reweighter

    def __call__(self, x):
        x = list(x)
        ref_distribution, reweighted_distribution = \
            self.reweighter.get_distributions(self.sim_data, self.ref_data, x)
        return np.sum(
            np.array(ref_distribution - reweighted_distribution) ** 2
        )


class DihedralOptimizer:
    """
    An CHARMM dihedral optimizer that allows you to change the member
    (depending on the reweighter provided) of multiplicities in the original FF
    """
    def __init__(
        self, sim_dihedrals, ref_dihedrals, reweighter,
        method='BFGS', options={'eps': 1e-2, 'gtol': 1e-04}
    ):
        """
        Args:
            sim_dihedrals: array, the simulation dihedral data (for each,
            not the average)
            ref_dihedrals: array, the reference state dihedral data
            reweighter: the CharmmDihedralReweighter
            method: string, the scipy.optimize routine
            options: dict, the options for the scipy optimizer
        """
        self.obj_func = ObjfuncDihedral(
            sim_dihedrals, ref_dihedrals, reweighter
        )
        self.reweighter = reweighter
        self.sim_dihedrals = sim_dihedrals
        self.ref_dihedrals = ref_dihedrals
        self.method = method
        self.options = options

    @property
    def num_multiplicity(self):
        return self.reweighter.multiplicity

    def __call__(self):
        # the dimension of the optimization is bound to the reweighter,
        # which is provided to the objective function
        start_k = np.array(self.reweighter.dihfunc2.k)
        optimum = sopt.minimize(
            self.obj_func, start_k, method=self.method, options=self.options
        )
        return optimum


def do_torsion_optmization(
    atoms, psf_file, crd_file, parameter_files,
    traj_to_compare, first_comp, last_comp,
    traj_to_fix, first_fix, last_fix,
    last_torfix=0,
    allowed_m=[1,2,3,4,5,6],
    temperature=323.15, nbins=100, plot=True
):
    """

    Args:
        atoms: list of strings, the four atom names in the dihedral
        psf_file: string, the psf file
        crd_file: string, the crd file
        parameter_files: list of the parameter files
        traj_to_compare: string, trajectory template of the reference state
        first_comp: integer, first trajectory of the reference state
        last_comp: integer, last trajectory of the regerence state
        traj_to_fix: string, trajectory template of the simulation
        first_fix: integer, first trajectory of the simulation the user want
        to use for the optimization
        last_fix: integer, last trajectory of the simulation the user
        want to use for the optimization
        allowed_m: list of integers, the multiplicity allowed in the
        optimization
        temperature: float, temperature of the simulation
        nbins: integer, number of bins for distribution of the dihedral
        plot: bool, if plot the distribution or not
    Returns:
        A dictionary use the multiplicity as key and
        [phase, force_constant] as content
    """
    from matplotlib import pyplot as plt
    target = DihedralTarget(atoms, psf_file, crd_file, parameter_files,
                             torsionfix=last_torfix)
    target.get_cosine_series()
    dihdata_fix = target.get_dihedrals(traj_to_fix, first_fix, last_fix)
    dihdata_ref = target.get_dihedrals(traj_to_compare, first_comp, last_comp)
    dih_f_old = dihedral_function(target.k, target.multp, target.phase)
    existing_ks = copy.deepcopy(dih_f_old.k)
    existing_ms = copy.deepcopy(dih_f_old.m)
    existing_ps = copy.deepcopy(dih_f_old.p)
    existing_ks_dict = {}
    existing_ps_dict = {}
    print("Here are the torsional terms before fitting:")
    for ek, em, ep in zip(existing_ks, existing_ms, existing_ps):
        print(
            str(ek) + '(kj/mole) OR ' +
            str(round(ek/4.184, 4)) +
            '(kcal/mole)', ', ', em, ', ', ep
        )
        existing_ks_dict[em] = ek
        existing_ps_dict[em] = ep
    new_ks = list(np.zeros(6))
    new_ms = list(np.arange(1, 7))
    new_ps = list(np.zeros(6))
    # Generate the new k and m lists that include more multiplicities
    for m in allowed_m:
        if m in existing_ms:
            new_k = existing_ks_dict[m]
            new_p = existing_ps_dict[m]
        else:
            # if missing the multiplicity,
            # use 0 for the force constant at beginning
            new_k = 0.0
            new_p = 0.0
        new_ks[m - 1] = new_k
        new_ps[m - 1] = new_p
        dih_f_free = dihedral_function(new_ks, new_ms, new_ps)
    # first try the fixed-m optimizer
    reweighter = CharmmDihedralReweighter(
        dih_f_old, dih_f_old, temperature=temperature
    )
    optimizer = DihedralOptimizer(dihdata_fix, dihdata_ref, reweighter)
    optim = optimizer()
    print("The residue from the fixed multiplicity fitting is:", optim.fun)
    # The following condition should be changed to judge the optim found above
    # is good or not
    two_stages = False
    if optim.fun > 0.00001 * nbins:
        two_stages = True
        print("Using more multiplicities ...")
        reweighter2 = CharmmDihedralReweighter(
            dih_f_old, dih_f_free, temperature=temperature
        )
        optimizer2 = DihedralOptimizer(dihdata_fix, dihdata_ref, reweighter2)
        optim2 = optimizer2()
    # prepare return value before plotting
    return_dic = {}
    print("Here are the fixed torsional terms:")
    if not two_stages:
        for tm, tp, tk in zip(target.multp, target.phase, list(optim.x)):
            print(
                str(tk) + '(kj/mole) OR ' +
                str(round(tk/4.184, 4)) +
                '(kcal/mole)', ', ', tm, ', ', tp
            )
            return_dic[tm] = [tp, tk]
    else:
        print('Before adding more ')
        for tm, tp, tk in zip(new_ms, new_ps, list(optim2.x)):
            print(
                str(tk) + '(kj/mole) OR ' + str(round(tk/4.184, 4)) +
                '(kcal/mole)', ', ', tm, ', ', tp
            )
            return_dic[tm] = [tp, tk]
    if plot:
        print('\n' + 'Plotting...')
        obj_func = ObjfuncDihedral(dihdata_fix, dihdata_ref, reweighter)
        ref_distrib, fixed_distrib = obj_func.reweighter.get_distributions(
            dihdata_fix, dihdata_ref, optim.x
        )
        _, unfixed_distrib = obj_func.reweighter.get_distributions(
            dihdata_fix, dihdata_ref, target.k
        )
        if two_stages:
            obj_func2 = ObjfuncDihedral(
                dihdata_fix, dihdata_ref, reweighter2
            )
            _, fixed_distrib2 = obj_func2.reweighter.get_distributions(
                dihdata_fix, dihdata_ref, optim2.x
            )
        xaxis = np.arange(-180, 180, 3.6)
        plt.figure(figsize=(12,7))
        plt.plot(
            xaxis, ref_distrib, label='Reference (PME)',
            alpha=0.75, linewidth=5, color='red'
        )
        plt.plot(
            xaxis, unfixed_distrib, label='Before fitting',
            alpha=0.75, linewidth=3, color='royalblue'
        )
        plt.plot(
            xaxis, fixed_distrib, label='After fitting',
            alpha=0.75, linewidth=5, color='darkorange'
        )
        if two_stages:
            plt.plot(
                xaxis, fixed_distrib2, label='Adding multiplicities',
                alpha=0.75, linewidth=5, color='yellowgreen'
            )
        plt.title(
            '{}-{}-{}-{}'.format(
                atoms[0].upper(), atoms[1].upper(),
                atoms[2].upper(), atoms[3].upper()
            ),
            fontsize = 24, fontname = 'URW Gothic'
        )
        plt.xticks(fontsize = 22, fontname = 'URW Gothic')
        plt.yticks(fontsize = 22, fontname = 'URW Gothic')
        plt.xlabel('Dihedral Angle', fontsize = 24, fontname = 'URW Gothic')
        plt.ylabel('Population', fontsize = 24, fontname = 'URW Gothic')
        plt.legend(fontsize = 18)
        plt.savefig(
            '/u/alanyu/jplots/tfix-{}-{}-{}-{}.png'.format(
            atoms[0], atoms[1], atoms[2], atoms[3]),
            dpi = 200, bbox_inches = 'tight'
        )
        plt.show()
    return return_dic



