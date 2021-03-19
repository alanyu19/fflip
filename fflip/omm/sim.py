# -*- coding: utf-8 -*-

"""Classes for Gromacs calculations"""

from __future__ import absolute_import, division, print_function

import shutil

from fflip.core.globconf import CONFIG

import fflip.core.fflipdir
from fflip.core.decorators import args_from_configfile
from fflip.core import shell
from fflip.omm import util as ommutil
from fflip.omm import exceptions as ommexception

class OmmCalculation(fflip.core.fflipdir.CoffeWorkDir):
    """
    Class for Openmm Charmm calculations (simulations//property/energy
    calculations//more_to_come...)
    """
    @fflip.core.fflipdir.log_exceptions
    @args_from_configfile
    def __init__(self, structure, psf_file, options_file, work_dir=".",
                 options={}, overwrite=False, checkpoint=None):
        """Supports :func:`~fflip.core.decorators.args_from_configfile`.
        Args:
            structure: structure file (.crd, .pdb,...). it does not have to
            exist upon construction.
            topology: psf file(.psf)
            work_dir: working directory (default: current working directory)
            overwrite: a boolean that specifies if a calculation is to be rerun,
            even if it was already terminated

        Raises:
            AssertionError: if the input arguments do not match
            OmmSimError: if something goes wrong in OpenMM
        """

        if overwrite:
            try:
                shutil.rmtree(work_dir)
            except OSError:
                pass

        super(OmmCalculation, self).__init__(
            work_dir, "OmmCalculation", locals()
        )

        if overwrite:
            self.logger.info(
                "Runs in overwrite mode: Old directory contents were removed."
            )
            # TODO: missing removing action here
        self.psf_file = self.abspath(psf_file, check_exists=True)
        self.structure = self.abspath(structure, check_exists=True)
        # DO NOT remove the following two lines!
        options['psf'] = psf_file
        options['crd'] = structure
        self.options = options
        self.options_file = self.abspath(options_file, check_exists=False)

        #self.checkpoint = None
        #if checkpoint is not None:
        #    self.checkpoint = self.abspath(checkpoint, check_exists=False)
        #self.terminated_file = os.path.join(self.fflip_dir, "terminated.txt")

        if self.options != {}:
            ommutil.gen_omm_options_file(self.options, self.options_file)


    @fflip.core.fflipdir.log_exceptions
    def __call__(self, command):
        """
        Run OpenMM Calculation.
        Args:
            command: str, the command to run the calculation
        Raises:
            OmmCalcError: If something goes wrong in OpenMM (to do)
        """
        print(self.logger)
        self.logger.info("Running calculation")
        self._ommrun(command)
        self.logger.info("Omm calculation finished.")

    @fflip.core.fflipdir.log_exceptions
    def _ommrun(self, command):
        """Run OpenMM CHARMM simulation from Rickflow."""
        try:
            self.call_cmd(command)
        except shell.ShellError as e:
            raise ommexception.OmmSimError(e, self.last_errfile)
