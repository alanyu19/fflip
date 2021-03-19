#!/usr/bin/env python

# -*- coding: utf-8 -*-


import os
import sys
from coffe.omm.util import get_md_options as gmd
from fflip.chm import *

# this can be passed as input but not neccessary for now
option_file_name = 'obscalc.inp'

opts = gmd(option_file_name)

first_dcd = int(sys.argv[1])

phos_pairs = ["Os1-HW", "O2-HW", "Os1-OW", "O2-OW"]

es_pairs = ["Ob-OW", "Os2-OW", "Ob-HW", "Os2-HW"]

os.system(
    "sbatch --array=0-{} calc.py {}".format(
        len(es_pairs + phos_pairs),
        first_dcd
    )
)
# 0 for the original energy
