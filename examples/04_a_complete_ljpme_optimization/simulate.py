#!/usr/bin/env python

from targets import *

root = '/u/alanyu/c36ljpme/fflow/'

sfile = '/u/alanyu/c36ljpme/fflow/le_no_o2l/solutions/gg/gg3/solution.txt'
tfile = '/u/alanyu/c36ljpme/fflow/le_no_o2l/solutions/gg/gg3/torfix.py'

sim_loc = '/v/gscratch/mbs/alanyu/c36ljpme'

# assert os.path.isfile(solfile)

index = range(39)

for i in index:
    properties[i].simulate(
        iteration='safe', trj_folder=sim_loc,
        change_para=True, torfix_file=tfile, solution_file=sfile
    )
