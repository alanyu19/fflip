#!/usr/bin/env python
# coding: utf-8

import os
import time
import numpy as np
import pandas as pd


class FutureResult(object):
    @staticmethod
    def get_result(
            path_to_file, file_name, time_to_wait, loading_method, delete=False
    ):
        while not os.path.isfile(path_to_file + file_name):
            time.sleep(time_to_wait)
        if loading_method == 'numpy':
            result = np.genfromtxt(path_to_file + file_name)
        elif loading_method == 'pandas':
            result = pd.read_csv(path_to_file + file_name)
        else:
            with open(path_to_file + file_name, 'r') as file:
                result = file.readlines()
        if delete:
            os.system("rm -rf {}".format(path_to_file))
        return result


def on_cluster(executable, executable_args_list, *args, **kwargs):
    out_dir = kwargs['out_dir']
    if os.path.isdir(out_dir):
        os.system("rm -rf {}".format(out_dir))
    os.system("mkdir {}".format(out_dir))
    if 'job_time' in kwargs:
        run_time = kwargs['job_time']
    else:
        run_time = '01:00:00'
    if 'partition' in kwargs:
        partition = kwargs['partition']
    else:
        partition = "ivy,sbr,hpodall,spodall"
    if 'ntasks' in kwargs:
        ntasks = kwargs['ntasks']
    else:
        ntasks = 1
    if 'conda_env' in kwargs:
        conda = kwargs['conda_env']
    else:
        conda = 'drude'
    with open(kwargs["submit_script"], 'w+') as f:
        f.write(
            "#!/bin/bash\n" +
            "#SBATCH --output=./{}/{}.out\n".format(
                out_dir, kwargs["slurm_name"]
            ) +
            "#SBATCH --error=./{}/{}.err\n".format(
                out_dir, kwargs["slurm_name"]
            ) +
            "#SBATCH --time={}\n".format(run_time) +
            "#SBATCH --partition={}\n".format(partition) +
            "#SBATCH --ntasks={}\n".format(ntasks) + '\n'
        )
        f.write(
            "source ~/.bashrc\n" +
            "conda activate {}\n".format(conda)
        )
        argstring = ''
        for arg in executable_args_list:
            argstring += ' {}'.format(str(arg))
        argstring += ' {}'.format(kwargs['out_dir'])
        f.write("\n" + "python " + executable + " " + argstring +
                " >& ./{}/{}.out".format(out_dir, kwargs["exec_name"]))
    os.system('sbatch {}'.format(kwargs['submit_script']))
    time.sleep(1)
    os.system("rm -f {}".format(kwargs['submit_script']))


