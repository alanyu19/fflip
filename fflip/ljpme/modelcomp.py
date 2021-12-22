# -*- coding: utf-8 -*-

# Contains master class of training target (area, rdf, scd, db ...)

import glob
import mdtraj as md
from openmm.app import NoCutoff, CharmmCrdFile, CharmmPsfFile, \
    CharmmParameterSet, HBonds
from openmm import Context, LangevinIntegrator
import openmm.unit as u

from fflip.omm.genclac import OmmJobGenerator
from fflip.omm.torsionfuncs import *
from fflip.chm import *
from fflip.drude import *
from fflip.omm.util import *


class ModelCompoundPool(object):
    def __init__(
        self, name, dihedrals, psf_file, crd_folder, box, temperature,
        last_seqno, sim_template=None, integrator=None, ff='additive',
        overwrite=True, start=True
    ):
        self.name = name
        self.dihedrals = dihedrals
        self.psf_file = psf_file
        self.crd_folder = crd_folder
        self.box = box
        self.temperature = temperature
        self.last_seqno = last_seqno
        self.sim_template = sim_template
        self.integrator = integrator
        self.ff = ff
        self.overwrite = overwrite

    def generate_model_compound(self):
        self.mc = ModelCompound(
            self.name, self.psf_file, self.dihedrals,
            self.sim_template, self.ff
        )

    def simulate(self, trj_folder, toppar_path, start):
        crd_files = glob.glob(os.path.join(self.crd_folder, '*'))
        for crd in crd_files:
            folder = os.path.join(trj_folder, crd.split('/')[-1])
            self.mc.simulate(
                folder, crd, toppar_path, self.box, self.temperature, self.last_seqno,
                self.integrator, self.overwrite, start
            )

    def dihedral_distribution(self, trj_folder):
        psf = CharmmPsfFile(self.psf_file)
        topology = md.Topology.from_openmm(psf.topology)
        trajs = []
        mc = self.mc
        # TODO: loop to include more seqno
        replicas = glob.glob(f"{trj_folder}/*")
        for rep in replicas:
            traj_file = os.path.join(
                rep, 'trj/dyn2.dcd'
            )
            trajs.append(md.load_dcd(traj_file, top=self.psf_file))
        for dihe in self.dihedrals:
            for i, traj in enumerate(trajs):
                analyzer = DihedralAngle(
                    topology, dihe[0], dihe[1], dihe[2], dihe[3]
                )
                a = analyzer(traj)
                dens, bins = np.histogram(
                    a, bins=360, range=(-np.pi, np.pi), density=True
                )
                if i == 0:
                    dens_all = dens
                else:
                    dens_all = dens + dens_all
            dens_all = dens_all / i
            to_save = np.array(
                [bins[1:] - bins[1] / 2 + bins[0] / 2, dens_all]
            ).transpose()
            np.savetxt(trj_folder + '/{}-{}-{}-{}.dat'.format(
                dihe[0], dihe[1], dihe[2], dihe[3]
            ), to_save)

        
class ModelCompound(Lipid, DrudeLipid):
    def __init__(
        self, name, psf_file, dihedrals=None, sim_template=None, ff='additive',
        charge_groups=None, lj_groups=None, nbtholes=None, charmm_group_list=[]
    ):
        self.name = name
        self.dihedrals = dihedrals
        self.psf_file = psf_file
        # self.crd_file = crd_file
        self.sim_template = sim_template
        self.ff = ff
        self.psf = CharmmPsfFile(self.psf_file)
        self.parameters = None
        self.topology = md.Topology.from_openmm(self.psf.topology)
        if self.ff is 'drude':
            DrudeLipid.__init__(self, name, charge_groups, lj_groups, nbtholes)
        elif self.ff is 'additive':
            Lipid.__init__(self, name=name, charmm_group_list=charmm_group_list)

    def simulate(self, trj_folder, crd, toppar_path, box, temperature, last_seqno,
                 integrator=None, overwrite=False, start=False, verbose=0):
        if self.sim_template is None:
            return 0  # exit

        trj_loc = os.path.realpath(trj_folder)
        if verbose >= 1:
            print("Runnning in {}".format(trj_loc))
        calc = OmmJobGenerator(
            crd, self.psf_file,
            template=self.sim_template, work_dir=trj_loc
        )
        # constructing options dictionary
        options = dict()
        # special treatment to the option file name
        options["mdo"] = "input.mdo"
        options["toppar_path"] = toppar_path
        if box is not None:
            options["box"] = float(box)
        else:
            warnings.warn("Better provide a box size, using 50 Å by default")
            options["box"] = 50
        options["temperature"] = temperature
        if integrator is not None:
            options["integrator"] = integrator
        else:
            options["integrator"] = "L"  # Langevin
        options["psf"] = self.psf_file
        options["crd"] = crd
        options["molecule"] = self.name
        options["ff"] = self.ff
        if last_seqno is not None:
            options["last_seqno"] = int(last_seqno)
        else:
            options["last_seqno"] = 2
        # get the job run
        job = calc(0, options, overwrite=overwrite)
        if start:
            job("rflow submit sdyn.sh")
    
    def load_parameter_files(self, parameter_files):
        """
        parameter_files: list of CHARMM parameter files
        """
        self.parameters = CharmmParameterSet(*parameter_files)
    
    def load_nbthole(self, params, offsets):
        for p, o in zip(params, offsets):
            if p.par_type is 'nbthole':
                self.add_nbthole(
                    p.center_names[0], p.center_names[1], p.original_p + o
                )
                
    def add_nbthole(self, atype1, atype2, nbt_value):
        self.parameters.nbthole_types[(atype1, atype2)] = nbt_value
        self.parameters.atom_types_str[atype1].nbthole[atype2] = nbt_value
        self.parameters.atom_types_str[atype2].nbthole[atype1] = nbt_value
        
    def generate_gas_system(self, parameter_files=None):
        if not self.parameters:
            assert parameter_files is not None
            self.load_parameter_files(parameter_files)
        self.system = self.psf.createSystem(
            self.parameters, 
            nonbondedMethod=NoCutoff,
            constraints=HBonds
        )
        
    def change_nb_params(self, params, offsets):
        for p, o in zip(params, offsets):
            if p.par_type is not 'nbthole':
                change_drude_ff_parameters(
                    self.system, self.topology, p, o, self.psf
                )
    
    def generate_context(self):
        self.context = Context(self.system, LangevinIntegrator(1, 1, 1))
        
    def generate_energies(self, crd_files):
        energies = []
        for crd_file in crd_files:
            crd = CharmmCrdFile(crd_file)
            self.context.setPositions(crd.positions)
            e = self.context.getState(getEnergy=True).getPotentialEnergy(
            ).value_in_unit(u.kilocalories_per_mole)
            energies.append(e)
        return np.array(energies)

