# -*- coding: utf-8 -*-

import os
def check_and_make_dir(directory: object) -> object:
    if not os.path.isdir(directory):  # directory for slurm output files
        os.mkdir(directory)

