#!/bin/bash
#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --partition=ivy,sbr,spodall,hwell,hpodall
#SBATCH --ntasks=1

# $1 is the path of trajectory files
# $SLURM_ARRAY_TASK_ID is the traj index 
RES=$(sbatch calc.sh $1 $SLURM_ARRAY_TASK_ID $2) && echo ${RES##* } >> jobs.id
