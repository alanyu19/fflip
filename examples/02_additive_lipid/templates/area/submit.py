# -*- coding: utf-8 -*-

#!/usr/bin/env python

import os
from fflip.omm.util import get_md_options as gmd


if not os.path.isdir('log'):
    os.mkdir('log')

if not os.path.isdir('block_data'):
    os.mkdir('block_data')

# this can be passed as input but not neccessary for now
option_file_name = 'obscalc.inp'

opts = gmd(option_file_name)

first_dcd = int(opts['first_trj'])
last_dcd = int(opts['last_trj'])
blk_size = int(opts['block_size']) # number of traj files each block

array_str = ""
for index, t in enumerate(range(first_dcd, last_dcd, blk_size)):
    if index == 0:
        array_str += str(t)
    else:
        array_str += ",{}".format(str(t))
print(array_str)
os.system("sbatch --array={} calc.py".format(array_str))
