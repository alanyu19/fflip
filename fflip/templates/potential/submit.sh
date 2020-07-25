#!/bin/bash
#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=sbr,sky,hpodall,ivy,spodall
#SBATCH --ntasks=1

# $1 is the path of trajectory files
# $SLURM_ARRAY_TASK_ID is used for selecting the target parameters
cwd=${PWD##*/}
RES=$(sbatch --job-name=$cwd calc.sh $1 $SLURM_ARRAY_TASK_ID $2 $3) && echo ${RES##* } >> jobs.id
