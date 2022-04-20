#!/bin/bash
##SBATCH --job-name *******************
#SBATCH -p a100
#SBATCH -t 12:00:00
#SBATCH --ntasks=1 --nodes=1 --qos=short
#SBATCH --gres=gpu:1

cd $SLURM_SUBMIT_DIR

module load cuda/11.0.3

cd $SLURM_SUBMIT_DIR

sleep 10

# run simulation and resubmit script
# module load cuDNN/7.0/CUDA-9.0

./dyn.py && rflow submit sdyn.sh && sleep 60
