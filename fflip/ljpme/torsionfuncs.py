# -*- coding: utf-8 -*-


from simtk.openmm import PeriodicTorsionForce
from simtk.openmm.app import LJPME

from rflow.trajectory import *
from fflip.omm.util import *


class DihedralFunction(object):
    def __init__(self, k, m, p):
        """
        Args:
            k: the force constant (in kcal / mole)
            m: multiplicity (integer, 1 ~ 6)
            p: phase (in rad)
        """
        self._k = k
        self._m = m
        self._p = p
        
    @property
    def k(self):
        return self._k
    
    @property
    def m(self):
        return self._m
    
    @property
    def p(self):
        return self._p
    
    def __call__(self, phi, input_in_rad=True):
        """
        Returns the energy for a given dihedral phi
        """
        if not input_in_rad:
            phi = np.pi * phi / 180
        e = 0
        for forc, mult, phas in zip(self.k, self.m, self.p):
            e += forc * (1 + np.cos(mult * phi - phas))
        return e


class CharmmDihedralReweighter(object):
    def __init__(self, dihfunc1, dihfunc2, temperature):
        """
        dihedral_funtion: the dihedral function object
        """
        self.dihfunc1 = dihfunc1
        self.dihfunc2 = dihfunc2
        self.temp = temperature
        
    @property
    def num_original_multiplicity(self):
        return len(list(self.dihfunc1.m))
    
    def calc_sim_energy(self, dihdata):
        return self.dihfunc1(dihdata)
    
    def sensitivity_per_k(self, dihdata, perturbation_of_k=0.02):
        """
        Calculate the sensitvity for each of the force constant,
        it's not finished and used but we might return to this
        at some point --YL
        Args:
            dihdata: series of dihedral angles, can be list or array
            perturbation_of_k: default=0.02, the perturbation of k used in
            the sensitivity calculation.
        Returns:
        """
        energy_list = []
        existing_ks = copy.deepcopy(self.dihfunc1.k)
        existing_ms = copy.deepcopy(self.dihfunc1.m)
        existing_ps = copy.deepcopy(self.dihfunc1.p)
        existing_ks_dict = {}
        existing_ps_dict = {}
        for ek, em, ep in zip(existing_ks, existing_ms, existing_ps):
            existing_ks_dict[em] = ek
            existing_ps_dict[em] = ep
        new_ks = list(np.zeros(6))
        new_ms = list(np.arange(1, 7))
        new_ps = list(np.zeros(6))
        # Generate the new k and m lists that cover all the 6 multiplicities
        for m in range(1, 7):
            if m in existing_ms:
                new_k = existing_ks_dict[m] + perturbation_of_k
                new_p = existing_ps_dict[m]
            else:
                new_k = perturbation_of_k
                new_p = 0.0
            new_ks[m-1] = new_k
            new_ps[m-1] = new_p
            perturbed_dih_func = DihedralFunction(new_ks, new_ms, new_ps)
            pe = perturbed_dih_func(dihdata)
            energy_list.append(pe)
        return energy_list
    
    def recalc_energy(self, data, k):
        new_dihfunc = DihedralFunction(k, self.dihfunc2.m, self.dihfunc2.p)
        return new_dihfunc(data)
        
    def get_distributions(self, current_ensemble_data, ref_ensemble_data, k):
        # get the energy for sim
        e_sim = self.calc_sim_energy(current_ensemble_data)
        # recalc the perturbed energy
        e_new = self.recalc_energy(current_ensemble_data, k)
        distrib_a = prob_distribution(
            ref_ensemble_data,
            100, -np.pi, np.pi,
            temperature=self.temp,
            original_energy=None,
            perturbed_energy=None,
            do_reweighting=False
        )
        distrib_b = prob_distribution(
            current_ensemble_data,
            100, -np.pi, np.pi,
            temperature=self.temp,
            original_energy=e_sim,
            perturbed_energy=e_new,
            do_reweighting=True
        )
        return distrib_a, distrib_b


class DihedralAngle(object):
    """
    Calculate the dihedral (average over all molecules) of four consecutive
    atoms
    """
    def __init__(self, topology, atom1, atom2, atom3, atom4):
        self.topology = topology
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.sele1 = self.topology.select("name {}".format(self.atom1))
        self.sele2 = self.topology.select("name {}".format(self.atom2))
        self.sele3 = self.topology.select("name {}".format(self.atom3))
        self.sele4 = self.topology.select("name {}".format(self.atom4))
        
    def __call__(self, traj):
        """
        Requirement:
        mdtraj installed and imported as md
        """
        dihedral_every_frame = md.compute_dihedrals(
            traj,
            np.swapaxes(
                np.array([self.sele1, self.sele2, self.sele3, self.sele4]), 0, 1
            )
        )
        return dihedral_every_frame
       

class DihedralTarget(object):
    def __init__(
        self, atoms, psf_file, parameter_files, torsionfix=0
    ):
        """
        the dihedral to calculate energy for and do reweighting on
        Args:
            atoms: list of the 4 atoms in the dihedral
            psf_file: psf file
            parameter_files: list of parameter files
            torsionfix: if the simulation is done on original + some_torsion_fix,
            then this torsion fix is needed for calculating the original energy
        """
        self.atoms = atoms
        self.psf_file = psf_file
        self.parameter_files = parameter_files
        self.torsionfix = torsionfix
    
    @property
    def name(self):
        for i, a in enumerate(self.atoms):
            if i == 0:
                name = a,
            else:
                name = name + '-' + a
        return name
    
    @property
    def psf(self): 
        psf, topology, parameters = read_structure_parameter_files(
            self.psf_file, self.parameter_files
        )
        return psf
    
    @property
    def top(self):
        psf, topology, parameters = read_structure_parameter_files(
            self.psf_file, self.parameter_files
        )
        return topology

    def create_system(self):
        """
        Create the system containing the torsion(s)
        Returns:
            the system containing the torsion(s)
        """
        cps = CharmmParameterSet(*self.parameter_files)
        x_y_length = 100 * u.angstrom
        z_length = 100 * u.angstrom
        psf = self.psf
        psf.setBox(x_y_length, x_y_length, z_length)
        system = psf.createSystem(
            cps,
            nonbondedMethod=LJPME,
            nonbondedCutoff=u.Quantity(value=1.0, unit=u.nanometer),
            constraints=HBonds,
            ewaldErrorTolerance=0.0005
        )
        if self.torsionfix != 0:
            # todo: test
            system = BrutalTorsionParameters(
                system, self.psf, self.torsionfix
            )
        self.system = system
        
    @property
    def torsion_force(self):
        if not hasattr(self, 'system'):
            self.create_system()
        for force in self.system.getForces():
            if isinstance(force, PeriodicTorsionForce):
                return force
    
    @property
    def _atomtypes(self):
        sele1 = self.top.select("name {}".format(self.atoms[0]))
        sele2 = self.top.select("name {}".format(self.atoms[1]))
        sele3 = self.top.select("name {}".format(self.atoms[2]))
        sele4 = self.top.select("name {}".format(self.atoms[3]))
        atoms = min((sele1[0], sele2[0], sele3[0], sele4[0]), 
                    (sele4[0], sele3[0], sele2[0], sele1[0]))
        for dih in self.psf.dihedral_list:
            _a1 = dih.atom1.idx
            _a2 = dih.atom2.idx
            _a3 = dih.atom3.idx
            _a4 = dih.atom4.idx
            _atoms = (_a1, _a2, _a3, _a4)
            if atoms == _atoms:
                break
        return dih.atom1.attype, dih.atom2.attype, dih.atom3.attype, \
            dih.atom4.attype

    def get_cosine_series(self):
        sele1 = self.top.select("name {}".format(self.atoms[0]))
        sele2 = self.top.select("name {}".format(self.atoms[1]))
        sele3 = self.top.select("name {}".format(self.atoms[2]))
        sele4 = self.top.select("name {}".format(self.atoms[3]))
        atoms = min((sele1[0], sele2[0], sele3[0], sele4[0]), 
                    (sele4[0], sele3[0], sele2[0], sele1[0]))
        count = 0
        k = []
        multp = []
        phase = []
        for torsion_id in range(self.torsion_force.getNumTorsions()):
            a1, a2, a3, a4, x, y, z = self.torsion_force.getTorsionParameters(
                torsion_id
            )
            _atoms = (a1, a2, a3, a4)
            if atoms == _atoms:
                count += 1
                if np.abs(z._value) > 0.00001:
                    multp.append(x)
                    phase.append(y._value)
                    k.append(z._value)  # in kJ/mole
        self.multp = multp
        self.phase = phase
        self.k = k
        self.k_kcal_per_mole = [kk * 0.239006 for kk in k]
        
    def get_atomtypes(self):
        """
        since the _atomtypes method/property takes several second to get,
        we want to recycle it (using 'atomtypes')
        """
        self.atomtypes = self._atomtypes
    
    def get_dihedrals(self, trj_template, first_trj, last_trj):
        trajs = TrajectoryIterator(
            first_sequence=first_trj, last_sequence=last_trj,
            filename_template=trj_template, 
            topology_file=self.psf_file,
            atom_selection="all", load_function=md.load_dcd
        )
        all_trj_data = []
        dihang = DihedralAngle(
            self.top,
            self.atoms[0], self.atoms[1], self.atoms[2], self.atoms[3]
        )
        for traj in trajs:
            data = dihang(traj)
            all_trj_data.append(data)
        return np.array(all_trj_data)


def prob_distribution(
    data, num_bins, lower_bound, upper_bound,
    temperature=323.15, original_energy=None, perturbed_energy=None,
    do_reweighting=False
):
    """
    data: can be any thing to get distribution
    num_bins: the number of bins we want
    lower_bound: the lower boundary of binning
    upper_bound: the upper boundary of binning
    temperature: the temeprature of the reweighting
    original_energy: in kJ/mole
    perturbed_energy: in kJ/mole
    do_reweighting:
    if set to False, return only the original distribution,
    if set to True, return both the original and the perturbed distributions.
    """
    data = np.array(data)
    width = upper_bound - lower_bound
    delta = width / num_bins
    bins = np.arange(lower_bound + delta, upper_bound + delta, delta)
    digitized = np.digitize(data, bins)
    original_distribution = [np.sum(digitized==i) for i in range(num_bins)]
    if not do_reweighting:
        return(original_distribution/np.sum(original_distribution))
    else:
        beta = beta_kjmol(temperature)
        o_energy = np.array(original_energy)
        p_energy = np.array(perturbed_energy)
        assert o_energy.shape == p_energy.shape == data.shape, \
        'Please make sure that the inputs (energies and the property data) ' \
        'are in the same shape! your data is {}, energies are {}(original) ' \
        'and {}(perturbed)'.format(data.shape, o_energy.shape, p_energy.shape)
        tune = beta * np.mean(p_energy) - beta * np.mean(o_energy)
        obs = [np.sum((digitized==i) * np.exp(
            -beta * p_energy + beta * o_energy + tune
        )) for i in range(num_bins)]
        ptf = np.sum(np.exp(-beta * p_energy + beta * o_energy + tune))
        return obs/ptf
    
    
def get_torsion_ids_slow(torsion_force, psf, at1, at2, at3, at4):
    id_list = []
    atoms = min((at1, at2, at3, at4), (at4, at3, at2, at1))
    for torsion_id in range(torsion_force.getNumTorsions()):
        _at1, _at2, _at3, _at4 = get_torsion_atomtypes(
            torsion_force, torsion_id, psf
        )
        _atoms = min((_at1, _at2, _at3, _at4), (_at4, _at3, _at2, _at1))
        if atoms == _atoms:
            id_list.append(torsion_id)
    return id_list


def gen_torsion_atomtypes_dic(torsion_force, psf):
    dihe_dic = {}
    for dih in psf.dihedral_list:
        a1=dih.atom1.idx
        a2=dih.atom2.idx
        a3=dih.atom3.idx
        a4=dih.atom4.idx
        assert a1 < a4
        dihe_dic[(a1, a2, a3, a4)] = (
            dih.atom1.type.name,
            dih.atom2.type.name,
            dih.atom3.type.name,
            dih.atom4.type.name
        )
    return dihe_dic


def get_torsion_ids_2(torsion_force, psf, at1, at2, at3, at4):
    id_list = []
    atoms = min((at1, at2, at3, at4), (at4, at3, at2, at1))
    dihe_dic = gen_torsion_atomtypes_dic(torsion_force, psf)
    for torsion_id in range(torsion_force.getNumTorsions()):
        a1, a2, a3, a4,_,_,_ = torsion_force.getTorsionParameters(torsion_id)
        assert a1 < a4
        _atoms = dihe_dic[(a1,a2,a3,a4)]
        if atoms == _atoms:
            id_list.append(torsion_id)
    return id_list


def get_torsion_ids(torsion_force, psf, at1, at2, at3, at4):
    atom_type_info = {}
    for a in psf.atom_list:
        atom_type_info[a.idx] = a.attype
    id_list = []
    atoms = min((at1, at2, at3, at4), (at4, at3, at2, at1))
    for torsion_id in range(torsion_force.getNumTorsions()):
        a1, a2, a3, a4,_,_,_ = torsion_force.getTorsionParameters(torsion_id)
        _atoms = min(
            (
                atom_type_info[a1], atom_type_info[a2],
                atom_type_info[a3], atom_type_info[a4]
            ),
            (
                atom_type_info[a4], atom_type_info[a3],
                atom_type_info[a2], atom_type_info[a1]
            )
        )
        if atoms == _atoms:
            id_list.append(torsion_id)
    return id_list


def get_torsion_atomtypes(torsion_force, torsion_id, psf):
    a1, a2, a3, a4, _, _, _ = torsion_force.getTorsionParameters(torsion_id)
    atoms = min((a1, a2, a3, a4), (a4, a3, a2, a1))
    for dih in psf.dihedral_list:
        _a1=dih.atom1.idx
        _a2=dih.atom2.idx
        _a3=dih.atom3.idx
        _a4=dih.atom4.idx
        _atoms = (_a1, _a2, _a3, _a4)
        if atoms == _atoms: #TODO: WILDCARDS
            break
    return dih.atom1.type.name, dih.atom2.type.name, dih.atom3.type.name, \
           dih.atom4.type.name


class TorsionFix(object):
    def __init__(self, atom1, atom2, atom3, atom4, values):
        self.atom1 = atom1.upper()
        self.atom2 = atom2.upper()
        self.atom3 = atom3.upper()
        self.atom4 = atom4.upper()
        self.values = values


def change_torsion(force, psf, tfs):
    for tf in tfs:
        idlist = get_torsion_ids(
            force, psf, tf.atom1, tf.atom2, tf.atom3, tf.atom4
        )
        oldmtps = {}
        for tid in idlist:
            idx1, idx2, idx3, idx4, mtp, phase, _ = \
            force.getTorsionParameters(tid)
            atmind = min([idx1, idx2, idx3, idx4],
                         [idx4, idx3, idx2, idx1])
            indstr = '{}-{}-{}-{}'.format(
                atmind[0], atmind[1], atmind[2], atmind[3]
            )
            if not indstr in oldmtps:
                oldmtps[indstr] = [mtp]
            else:
                oldmtps[indstr].append(mtp)
            if mtp in tf.values:
                new_torsion = [
                    idx1, idx2, idx3, idx4,
                    mtp, phase, u.Quantity(
                        tf.values[mtp][1], unit = u.kilojoule_per_mole
                    )
                ]
                force.setTorsionParameters(tid, *new_torsion)
        for indstr in oldmtps:
            for m in tf.values:
                if m not in oldmtps[indstr]:
                    indx = indstr.split('-')
                    force.addTorsion(
                        int(indx[0]), int(indx[1]), int(indx[2]), int(indx[3]),
                        m, tf.values[m][0], u.Quantity(
                        tf.values[m][1], unit = u.kilojoule_per_mole
                        )
                    )

            
def brutal_torsion_parameters(system, psf, tfixes):
    for force in system.getForces():
        if isinstance(force, PeriodicTorsionForce):
            break
    change_torsion(force, psf, tfixes)
    return system
