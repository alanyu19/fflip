#!/usr/bin/python

from fflip.analysis.edp_util import *

PSF_FILE = "../sim/plant.psf"
combine_to_average_edp(PSF_FILE)
combine_to_average(PSF_FILE)
