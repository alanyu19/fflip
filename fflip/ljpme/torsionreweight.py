# -*- coding: utf-8 -*-

import scipy.optimize as sopt
from fflip.ljpme.torsionfuncs import *
from fflip.ljpme.util import isnotebook


def prob_to_pmf_kcal_per_mole(dist, temperature, offset=0):
    dist = dist + offset  # offset can be used to avoid overflow when dist=0
    return (
        -np.log(dist/np.sum(dist)) - min(-np.log(dist/np.sum(dist)))
    ) * (6.02 * 1.3806 * temperature * 0.000239006)


class ObjfuncDihedral(object):
    def __init__(self, sim, ref, reweighter, temperature, criterion="prob",
                 method='ensemble', boltzmann=True, pmf_weight=1):
        """
        Args:
            sim (array like): the simulated dihedral distribution.
            ref (array like): the reference dihedral distribution.
            reweighter: a CharmmDihedralReweighter object.
            criterion (str): can be "pmf" or "prob" (probability).
            boltzmann: if use Boltzmann weighting for "pmf"
        """
        self.sim_data = sim
        self.ref_data = ref
        self.reweighter = reweighter
        self.temperature = temperature
        self.criterion = criterion
        self.method = method
        self.boltzmann = boltzmann
        self.pmf_weight = pmf_weight

    def __call__(self, x):
        x = list(x)
        ref_distribution, reweighted_distribution = \
            self.reweighter.get_distributions(
                self.sim_data, self.ref_data, x, method=self.method
            )
        if self.criterion == "prob":
            return np.sum(
                np.array(ref_distribution - reweighted_distribution)**2
            ) * ref_distribution.shape[0]  # to account the influence of num_bins
        elif self.criterion == "pmf" or self.criterion == "prob+pmf":
            ref_e_in_kcal_per_mole = prob_to_pmf_kcal_per_mole(
                ref_distribution, self.temperature, 0.00001
            )  # offset to deal with zero prob
            rew_e_in_kcal_per_mole = prob_to_pmf_kcal_per_mole(
                reweighted_distribution, self.temperature, 0.00001
            )  # offset to deal with zero prob
            if self.boltzmann:
                factor = self.pmf_weight * 0.5 * (ref_distribution + reweighted_distribution) * \
                    (ref_distribution != 0) * (reweighted_distribution != 0) * \
                    ref_distribution.shape[0]  # to account the influence of num_bins
            else:
                factor = self.pmf_weight * (ref_distribution != 0) * (reweighted_distribution != 0)
            ssr = np.mean(factor * (ref_e_in_kcal_per_mole - rew_e_in_kcal_per_mole)**2)
            if self.criterion == "prob+pmf":
                ssr += sum(
                    np.array(ref_distribution - reweighted_distribution)**2
                ) * ref_distribution.shape[0]  # to account the influence of num_bins
            return ssr
        else:
            raise Exception(
                "Criterion `{}` not accepted! Only support 'prob', 'pmf', or 'prob+pmf'".format(
                    self.criterion
                )
            )


class DihedralOptimizer(object):
    """
    An CHARMM dihedral optimizer that allows you to change the member
    (depending on the reweighter provided) of multiplicities in the original FF
    """
    def __init__(
        self, sim_dihedrals, ref_dihedrals, reweighter, temperature, criterion="prob",
        boltzmann=True, method='BFGS', options={'eps': 1e-2, 'gtol': 1e-04},
        pmf_weight=1, perturb_method='ensemble'
    ):
        """
        Args:
            sim_dihedrals: array, the simulation dihedral data (for each, not the average)
            ref_dihedrals: array, the reference state dihedral data
            reweighter: the CharmmDihedralReweighter
            temperature: float, temperature of the systems
            criterion: string, can be 'prob', 'pmf', and 'prob+pmf'
            boltzmann: bool, use boltzmann weighting in error function
            method: string, the scipy.optimize routine
            options: dict, the options for the scipy optimizer
            pmf_weight: float, the relative weight of pmf (compared to prob) when use 'prob+pmf'
            perturb_method: string, can be 'ensemble' and 'single'. If perturb the ensemble
            (exact but information can be washed out due to noise) or single molecules
            (approximation in the sense of stat-mech but more robust in practice)
        """
        self.obj_func = ObjfuncDihedral(
            sim_dihedrals, ref_dihedrals, reweighter, temperature, criterion=criterion,
            method=perturb_method, boltzmann=boltzmann, pmf_weight=pmf_weight
        )
        self.reweighter = reweighter
        self.temperature = temperature
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


def torsion_match_two_ff(
        atoms,
        psf_file_comp, parameter_files_comp,
        psf_file_fix, parameter_files_fix,
        traj_comp, first_comp, last_comp,
        traj_fix, first_fix, last_fix,
        last_torfix=0,
        allowed_m=[1, 2, 3, 4, 5, 6],
        criterion="prob", 
        boltzmann=True,
        temperature=323.15,
        min_ssr_to_add_mult=5,
        perturb_method='ensemble',
        pmf_weight=1,
        opt_method="BFGS",
        opt_options={'eps': 1e-2, 'gtol': 1e-04},
        plot=True, save_plot_to="."
):
    """
    Predict the good force constants to match the reference ensemble's dihedral
    distribution. Multiplicities are expanded to "allowed_m" if original multiplicities
    (read from the parameter files) can't fulfill the "min_ssr_to_add_mult".
    Args:
        atoms: list of strings, the four atom names in the dihedral
        psf_file_comp: string, the psf file of the comparison system
        parameter_files_comp: list of the parameter files
        psf_file_fix: string, the psf file of the system to do the dihedral matching
        traj_to_compare: string, trajectory template of the reference system
        first_comp: integer, first trajectory index of the reference system
        last_comp: integer, last trajectory index of the regerence system
        traj_fix: string, trajectory template of the system to the dihedral matching
        first_fix: integer, first trajectory index used for the optimization
        last_fix: integer, last trajectory index used for the optimization
        allowed_m: list of integers, the multiplicities allowed in the optimization
        criterion (str): can be "pmf" or "prob" (probability)
        boltzmann (bool): if use Boltzmann weighting for "pmf"
        temperature (float): temperature of the simulation
        nbins (integer): number of bins for distribution of the dihedral
        min_ssr_to_add_mult (float): minimum ssr to expand multiplicity
        plot (bool): if plot the distribution or not
    Returns:
        A dictionary use the multiplicity as key and
        [phase, force_constant] as content
    """
    from matplotlib import pyplot as plt
    from datetime import datetime
    # print('step1', datetime.now())
    target = DihedralTarget(
        atoms, psf_file_fix, parameter_files_fix, torsionfix=last_torfix
    )
    target_comp = DihedralTarget(
        atoms, psf_file_comp, parameter_files_comp, torsionfix=last_torfix
    )
    # print('step2', datetime.now())
    target.get_cosine_series()
    target_comp.get_cosine_series()
    # print('step3', datetime.now())
    dihdata_fix = target.get_dihedrals(traj_fix, first_fix, last_fix)
    dihdata_ref = target_comp.get_dihedrals(traj_comp, first_comp, last_comp)
    #  print('step4', datetime.now())
    dih_f_old = DihedralFunction(target.k, target.multp, target.phase)
    # print('step5', datetime.now())
    existing_ks = copy.deepcopy(dih_f_old.k)
    existing_ms = copy.deepcopy(dih_f_old.m)
    existing_ps = copy.deepcopy(dih_f_old.p)
    # print('finish', datetime.now())
    existing_ks_dict = {}
    existing_ps_dict = {}
    # print("Here are the torsional terms before fitting:\n")
    for ek, em, ep in zip(existing_ks, existing_ms, existing_ps):
        print(
            str(round(ek, 4)) + '(kj/mole) OR ' +
            str(round(ek/4.184, 4)) + '(kcal/mole)',
            ', for m={}'.format(em),
            ', phase={}'.format(round(180*ep/np.pi, 0))
        )
        existing_ks_dict[em] = ek
        existing_ps_dict[em] = ep
    print("\n\n")
    new_ks = list(np.zeros(len(allowed_m)))
    new_ms = allowed_m
    new_ps = list(np.zeros(len(allowed_m)))
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
        dih_f_free = DihedralFunction(new_ks, new_ms, new_ps)
    # first try the fixed-m optimizer
    reweighter = CharmmDihedralReweighter(
        dih_f_old, dih_f_old, temperature=temperature
    )
    optimizer = DihedralOptimizer(
        sim_dihedrals=dihdata_fix,
        ref_dihedrals=dihdata_ref,
        reweighter=reweighter,
        temperature=temperature,
        criterion=criterion,
        boltzmann=boltzmann,
        method=opt_method,
        options=opt_options,
        pmf_weight=pmf_weight,
        perturb_method=perturb_method
    )
    optim = optimizer()
    print("The residue from the fixed multiplicity fitting is:", optim.fun)
    # The following condition should be changed to judge the optim found above
    # is good or not
    two_stages = False
    if optim.fun > min_ssr_to_add_mult:
        two_stages = True
        print("Expanding multiplicities ...\n")
        reweighter2 = CharmmDihedralReweighter(
            dih_f_old, dih_f_free, temperature=temperature
        )
        optimizer2 = DihedralOptimizer(
            sim_dihedrals=dihdata_fix,
            ref_dihedrals=dihdata_ref,
            reweighter=reweighter2,
            temperature=temperature,
            criterion=criterion,
            boltzmann=boltzmann,
            pmf_weight=pmf_weight,
            perturb_method=perturb_method
        )
        optim2 = optimizer2()
    # prepare return value before plotting
    return_dic = {}
    print("Here are the fixed torsional terms:\n")
    if not two_stages:
        pass
    else:
        print('Before adding more multiplicities:')
    for tm, tp, tk in zip(target.multp, target.phase, list(optim.x)):
        print(
            str(round(tk, 4)) + '(kj/mole) OR ' +
            str(round(tk / 4.184, 4)) + '(kcal/mole)',
            ', for m={}'.format(tm),
            ', phase={}'.format(round(180 * tp / np.pi, 0)),
        )
        if not two_stages:
            return_dic[tm] = [tp, tk]
        else:
            pass
    print("\n")
    if two_stages:
        print("After adding more multiplicities:")
        for tm, tp, tk in zip(new_ms, new_ps, list(optim2.x)):
            print(
                str(round(tk, 4)) + '(kj/mole) OR ' +
                str(round(tk / 4.184, 4)) + '(kcal/mole)',
                ', for m={}'.format(tm),
                ', phase={}'.format(round(180 * tp / np.pi, 0))
            )
            return_dic[tm] = [tp, tk]
        print("\n\n")
    if plot:
        obj_func = ObjfuncDihedral(
            dihdata_fix, dihdata_ref, reweighter, temperature, method=perturb_method
        )
        ref_distrib, fixed_distrib = obj_func.reweighter.get_distributions(
            dihdata_fix, dihdata_ref, optim.x, method=perturb_method
        )
        _, unfixed_distrib = obj_func.reweighter.get_distributions(
            dihdata_fix, dihdata_ref, target.k, method=perturb_method
        )
        if two_stages:
            obj_func2 = ObjfuncDihedral(
                dihdata_fix, dihdata_ref, reweighter2, temperature, method=perturb_method
            )
            _, fixed_distrib2 = obj_func2.reweighter.get_distributions(
                dihdata_fix, dihdata_ref, optim2.x, method=perturb_method
            )
        xaxis = np.arange(-180, 180, 3.6)
        plt.figure(figsize=(6, 4))
        plt.plot(
            xaxis, ref_distrib, label='Reference',
            alpha=0.75, linewidth=4, color='red'
        )
        plt.plot(
            xaxis, unfixed_distrib, label='Before fitting',
            alpha=0.75, linewidth=2, color='k'
        )
        plt.plot(
            xaxis, fixed_distrib, label='After fitting',
            alpha=0.75, linewidth=2, color='b'
        )
        if two_stages:
            plt.plot(
                xaxis, fixed_distrib2, label='Adding multiplicities',
                alpha=0.75, linewidth=3, color='g'
            )
            plt.title(
            "{}-{}-{}-{}".format(
                atoms[0].upper(), atoms[1].upper(),
                atoms[2].upper(), atoms[3].upper()
            ), fontsize=18
        )
        plt.xticks(fontsize=14)  # fontname='URW Gothic'
        plt.yticks(fontsize=14)
        plt.xlabel('Dihedral Angle', fontsize=18)
        plt.ylabel('(Relative) Population', fontsize=18)
        plt.legend(fontsize=12)
        plt.savefig(
            os.path.join(save_plot_to, 'tfix-{}-{}-{}-{}.png'.format(
                atoms[0], atoms[1], atoms[2], atoms[3]
            )), dpi=200, bbox_inches='tight'
        )
        if isnotebook():
            plt.show()
        plt.close()
        
        # PMF plot
        plt.figure(figsize=(6, 4))
        # ref_e_in_kt = - np.log(ref_distribution + 0.0001) 
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plt.plot(
                xaxis, prob_to_pmf_kcal_per_mole(ref_distrib, temperature),
                label='Reference', alpha=0.75, linewidth=4, color='r'
            )
            plt.plot(
                xaxis, prob_to_pmf_kcal_per_mole(unfixed_distrib, temperature),
                label='Before fitting', alpha=0.75, linewidth=2, color='k'
            )
            plt.plot(
                xaxis, prob_to_pmf_kcal_per_mole(fixed_distrib, temperature),
                label='After fitting', alpha=0.75, linewidth=2, color='b'
            )
        if two_stages:
            plt.plot(
                xaxis, prob_to_pmf_kcal_per_mole(fixed_distrib2, temperature),
                label='Adding multiplicities', alpha=0.75, linewidth=3, color='g'
            )
            plt.title(
            "{}-{}-{}-{} (PMF)".format(
                atoms[0].upper(), atoms[1].upper(),
                atoms[2].upper(), atoms[3].upper()
            ), fontsize=18
        )
        plt.xticks(fontsize=14)  # fontname='URW Gothic'
        plt.yticks(fontsize=14)
        plt.xlabel('Dihedral Angle', fontsize=18)
        plt.ylabel('PMF (kcal/mole)', fontsize=18)
        plt.legend(fontsize=12)
        plt.savefig(
            os.path.join(save_plot_to, 'tfix-{}-{}-{}-{}-pmf.png'.format(
                atoms[0], atoms[1], atoms[2], atoms[3]
            )), dpi=200, bbox_inches='tight'
        )
        if isnotebook():
            plt.show()
        plt.close()
    return return_dic


def perturb_4_ctl2(psf_file, parameter_files,
                   traj_template, trj_index, sn2=(2, 15), sn1=(2, 15),
                   perturbation=0.01):
    print("Preparing Dihedral Targets ...")
    dt_list = []
    for i in range(sn2[0], sn2[1] - 2):
        dt = DihedralTarget(
            ['C2{}'.format(i), 'C2{}'.format(i+1),
             'C2{}'.format(i+2), 'C2{}'.format(i+3)],
            psf_file,
            parameter_files=parameter_files, torsionfix=0
        )
        if i==sn2[0]:
            dt.create_system()
            dt.get_cosine_series()
        else:
            dt.system = dt_list[0].system
            dt.multp = dt_list[0].multp
            dt.phase = dt_list[0].phase
            dt.k = dt_list[0].k
            dt.k_kcal_per_mole = dt_list[0].k_kcal_per_mole
        dt_list.append(dt)
    for i in range(sn1[0], sn1[1] - 2):
        dt = DihedralTarget(
            ['C3{}'.format(i), 'C3{}'.format(i+1),
             'C3{}'.format(i+2), 'C3{}'.format(i+3)],
            psf_file,
            parameter_files=parameter_files, torsionfix=0
        )
        dt.system = dt_list[0].system
        dt.multp = dt_list[0].multp
        dt.phase = dt_list[0].phase
        dt.k = dt_list[0].k
        dt.k_kcal_per_mole = dt_list[0].k_kcal_per_mole
        dt_list.append(dt)
    print('Finished')
    # start to get the dihedrals and energies
    tj = trj_index
    oelist = []
    pelist = [[] for _ in range(len(dt_list[0].k))]
    for dt in dt_list:
        if isinstance(trj_index, list):
            dd = dt.get_dihedrals(traj_template, tj[0], tj[1])
        elif isinstance(trj_index, int):
            dd = dt.get_dihedrals(traj_template, tj, tj)
        dih_f_old = DihedralFunction(dt.k, dt.multp, dt.phase)
        e_old = dih_f_old(dd)
        oelist.append(e_old)
        for ki in range(len(dt.k)):
            ptbd_k = list(np.array(dt.k))
            ptbd_k[ki] += perturbation
            dih_f_ptbd = DihedralFunction(ptbd_k, dt.multp, dt.phase)
            e_new = dih_f_ptbd(dd)
            pelist[ki].append(e_new)
    return oelist, pelist


def perturb_ctl2_ctl2_ctl2_ctl3(
    psf_file, parameter_files,
    traj_template, trj_index, sn2=16, sn1=16,
    perturbation=0.01
):
    print("Preparing Dihedral Targets ...")
    dt_list = []
    dt = DihedralTarget(
        ['C2{}'.format(sn2-3), 'C2{}'.format(sn2-2),
         'C2{}'.format(sn2-1), 'C2{}'.format(sn2)],
        psf_file,
        parameter_files=parameter_files, torsionfix=0
    )
    dt.create_system()
    dt.get_cosine_series()
    dt_list.append(dt)
    dt = DihedralTarget(
        ['C3{}'.format(sn1 - 3), 'C3{}'.format(sn1 - 2),
         'C3{}'.format(sn1 - 1), 'C3{}'.format(sn1)],
        psf_file,
        parameter_files=parameter_files, torsionfix=0
    )
    dt.create_system()
    dt.get_cosine_series()
    dt_list.append(dt)
    print('Finished')
    # start to get the dihedrals and energies
    tj = trj_index
    oelist = []
    pelist = [[] for _ in range(len(dt_list[0].k))]  # assume homogeneous
    for dt in dt_list:
        if isinstance(trj_index, list):
            dd = dt.get_dihedrals(traj_template, tj[0], tj[1])
        elif isinstance(trj_index, int):
            dd = dt.get_dihedrals(traj_template, tj, tj)
        dih_f_old = DihedralFunction(dt.k, dt.multp, dt.phase)
        e_old = dih_f_old(dd)
        oelist.append(e_old)
        for ki in range(len(dt.k)):
            ptbd_k = list(np.array(dt.k))
            ptbd_k[ki] += perturbation
            dih_f_ptbd = DihedralFunction(ptbd_k, dt.multp, dt.phase)
            e_new = dih_f_ptbd(dd)
            pelist[ki].append(e_new)
    return oelist, pelist


def perturb_any_dihedral(
    atoms, psf_file, parameter_files, torfix,
    traj_template, trj_index, perturbation=0.01):
    dt = DihedralTarget(
        atoms,
        psf_file,
        parameter_files=parameter_files,
        torsionfix=torfix
    )
    dt.create_system()
    dt.get_cosine_series()
    # start to get the dihedrals and energies
    tj = trj_index
    pe_dict = dict()
    if isinstance(trj_index, list):
        dd = dt.get_dihedrals(traj_template, tj[0], tj[1])
    elif isinstance(trj_index, int):
        dd = dt.get_dihedrals(traj_template, tj, tj)
    dih_f_old = DihedralFunction(dt.k, dt.multp, dt.phase)
    e_old = dih_f_old(dd)
    # oelist.append(e_old)
    for ki in range(len(dt.k)):
        ptbd_k = list(np.array(dt.k))
        mtps = list(np.array(dt.multp))
        ptbd_k[ki] += perturbation
        dih_f_ptbd = DihedralFunction(ptbd_k, dt.multp, dt.phase)
        e_new = dih_f_ptbd(dd)
        pe_dict[mtps[ki]] = e_new
    return e_old, pe_dict
