# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import os
from fflip.omm.sim import *
from fflip.omm.util import copy_folder, copy_one_layer 


class OmmJobGenerator:
    """A generator for omm/mdtraj simulation/calculation plans."""
    def __init__(self, structure, psf_file, template, work_dir):
        """
        Args:
            structure: str, crd/pdb
            psf_file: str, psf
            template: str, path of the template folder
            work_dir:
        """
        self.structure = structure
        self.psf_file = psf_file
        self.template = template
        self.work_dir = work_dir

    def gen_pot_calc(self, options, overwrite=False,
                     options_file='potcalc.inp'):
        """
        To generate a potential calculation based on the input parameters
        Args:
            options: dict, the input parameters.
            overwrite: bool, if overwrite the old result or not.
            options_file: str, a file stores the inputs.
        Returns:
            A portential calculation that can be called.

        """
        assert isinstance(options, dict)
        must_in = ['trj_location', 'first_trj', 'last_trj',
                   'block_size', 'option_file']
        for key in must_in:
            assert key in options, "{} must in options".format(key)
        job_template = self.template
        work_dir = self.work_dir
        if os.path.exists(work_dir):
            if overwrite == True:
                os.system('rm -r {}'.format(work_dir))
                copy_folder(job_template, work_dir)
            else:
                copy_one_layer(job_template, work_dir)
        else:
            copy_folder(job_template, work_dir)
        if 'percentage' not in options:
            # TODO: do we want to keep this percentage of parameter offset
            #  used in the sensitivity analysis in fflip?
            options['percentage'] = 1
        if 'solution' not in options:
            # TODO: same question here!
            options['solution'] = None
        elif options['solution'] == None:
            pass
        else:
            os.system(
                "cp {} {}/solution.txt".format(options['solution'], work_dir))
        if 'torfix' not in options:
            options['torfix'] = None
        elif options['torfix'] == None:
            pass
        else:
            os.system("cp {} {}/torfix.py".format(options['torfix'], work_dir))
        energy_calculation = OmmCalculation(
            self.structure, self.psf_file, options_file=options_file,
            work_dir=work_dir, options=options
        )
        return energy_calculation

    def gen_prop_calc(self, options, overwrite=False):
        """
        Generate Energy Calculation
        """
        assert isinstance(options, dict)
        must_in = ['trj_location', 'first_trj', 'last_trj',
                   'block_size', 'option_file']
        for key in must_in:
            assert key in options, "{} must in options".format(key)
        # copy over the template
        job_template = self.template
        work_dir = self.work_dir
        if os.path.exists(work_dir):
            if overwrite is True:
                os.system('rm -r {}'.format(work_dir))
                copy_folder(job_template, work_dir)
            else:
                copy_one_layer(job_template, work_dir)
        else:
            copy_folder(job_template, work_dir)
        property_calculation = OmmCalculation(
            self.structure, self.psf_file, options_file=options['option_file'],
            work_dir=work_dir, options=options
        )
        return property_calculation

    def gensim(self, options, overwrite=False):
        """
        To generate a simulation (OmmChmCalculation) based on
        some template and input parameters.
        Args:
            options: dict, the input parameters.
            template: str, the path of the template relative to
            the template root.
            overwrite: bool, if overwrite the old result or not.

        Returns:
            A portential calculation that can be called.
        """
        assert isinstance(options, dict)
        must_in = ['mdo', 'last_seqno']
        for key in must_in:
            assert key in options, "{} must in options".format(key)
        sim_template = self.template
        work_dir = self.work_dir
        if os.path.exists(work_dir):
            if overwrite == True:
                os.system('rm -r {}'.format(work_dir))
                copy_folder(sim_template, work_dir)
            else:
                copy_one_layer(sim_template, work_dir)
        else:
            copy_folder(sim_template, work_dir)
        with open(os.path.join(work_dir, 'last.seqno'), 'w') as f:
            f.write(str(options['last_seqno']))
        simulation = OmmCalculation(
            self.structure, self.psf_file,
            options_file=options['mdo'],
            work_dir=work_dir, options=options
        )
        # TODO: we added 'param_offset_file' and 'last_seqno' in options,
        #  check if these will influence our simulations ...
        return simulation

    def __call__(self, calc_type, options, overwrite=False):
        # TODO: iteration number should not appear in this function/class,
        #  but it should appear in the work directory, think about how to do
        #  this.
        if calc_type == 0: # simulation
            calculation = self.gensim(options, overwrite=overwrite)
            return calculation
        if calc_type == 1: # potential
            calculation = self.gen_pot_calc(options, overwrite=overwrite)
            return calculation
        if calc_type == 2: # observable
            calculation = self.gen_prop_calc(options, overwrite=overwrite)
            return calculation
