# -*- coding: utf-8 -*-

# Contains master class of training target (area, rdf, scd, db ...)

import glob
import mdtraj as md
from fflip.omm.genclac import OmmJobGenerator
from fflip.omm.torsionfuncs import *

class ModelCompoundPool(object):
    def __init__(
        self, name, dihedrals, psf_file, crd_folder, box, temperature, last_seqno, 
        sim_template=None, integrator=None, ff='additive', overwrite=True, start=True):
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
            self.name, self.dihedrals, self.psf_file,
            self.sim_template, self.ff
        )

    def simulate(self, trj_folder, start):
        crd_files = glob.glob(os.path.join(self.crd_folder, '*'))
        for crd in crd_files:
            folder = os.path.join(trj_folder, crd.split('/')[-1])
            self.mc.simulate(
                folder, crd, self.box, self.temperature, self.last_seqno,
                self.integrator, self.overwrite, start
            )

    def dihedral_distribution(self, trj_folder):
        psf = CharmmPsfFile(self.psf_file)
        topology = md.Topology.from_openmm(psf.topology)
        trajs = []
        for mc in self.mc_list:
            # TODO: loop to include more seqno
            traj_file = os.path.join(trj_folder, mc.crd_file.split('/')[-1], 'trj/dyn2.dcd')
            trajs.append(md.load_dcd(traj_file, top=self.psf_file))
        for dihe in self.dihedrals:
            for i, traj in enumerate(trajs):
                analyzer = DihedralAngle(topology, dihe[0], dihe[1], dihe[2], dihe[3])
                # ti = h[i]
                # index = f"{ti[0]+1}_{ti[1]+1}_{ti[2]+1}"
                a = analyzer(traj)
                dens, bins = np.histogram(a, bins=360, range=(-np.pi, np.pi), density=True)
                if i == 0:
                    dens_all = dens
                else:
                    dens_all = dens + dens_all
            dens_all = dens_all / i
            to_save = np.array([bins[1:] - bins[1] / 2 + bins[0] / 2, dens_all]).transpose()
            np.savetxt(trj_folder + '/{}-{}-{}-{}.dat'.format(dihe[0], dihe[1], dihe[2], dihe[3]), to_save)

        
class ModelCompound(object):
    def __init__(
        self, name, dihedrals, psf_file, sim_template=None, ff='additive'
    ):
        self.name = name
        self.dihedrals = dihedrals
        self.psf_file = psf_file
        self.sim_template = sim_template
        self.ff = ff

    def simulate(self, trj_folder, crd, box, temperature, last_seqno,
                 integrator=None,
                 overwrite=False, start=False, verbose=0):
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
        if box is not None:
            options["box"] = float(box)
        else:
            warnings.warn("Better provide a box size, using 50 Å by default")
            options["box"] = 50
        options["temperature"] = temperature
        if integrator is not None:
            options["intgrt"] = integrator
        else:
            options["intgrt"] = "L"  # Langevin
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
