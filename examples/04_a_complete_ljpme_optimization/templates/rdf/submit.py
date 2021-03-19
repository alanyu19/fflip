#!/usr/bin/env python

# -*- coding: utf-8 -*-


import os
from coffe.omm.util import get_md_options as gmd


if not os.path.isdir('log'):
    os.mkdir('log')

if not os.path.isdir('block_data'):
    os.mkdir('block_data')

# this can be passed as input but not neccessary for now
option_file_name = 'obscalc.inp'

opts = gmd(option_file_name)

first_dcd = int(opts['first_trj'])
last_dcd = int(opts['last_trj'])
blk_size = int(opts['block_size'])  # number of traj files each block

for index, t in enumerate(range(first_dcd, last_dcd + 1, blk_size)):
    os.system("python seedcalc.py {}".format(t))
