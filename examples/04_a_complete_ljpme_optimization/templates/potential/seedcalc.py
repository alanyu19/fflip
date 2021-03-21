# -*- coding: utf-8 -*-
#!/usr/bin/env python

import os
import sys
from coffe.omm.util import get_md_options as gmd
from fflip.chm import *

# this can be passed as input but not neccessary for now
option_file_name = 'potcalc.inp'

opts = gmd(option_file_name)

first_dcd = int(sys.argv[1])
blk_size = int(opts['block_size'])  # number of traj files each block
lipname = str(opts['lipname'])

lipid = parse_lipid(lipname)

parameter_sets = lipid.parse_gtcnp()

os.system(
    "sbatch --array=0-{} calc.py {}".format(len(parameter_sets), first_dcd)
)
# 0 for the original energy
