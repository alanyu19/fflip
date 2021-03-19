#!/bin/bash
#SBATCH --output=./log/master_%A_%a.out
#SBATCH --error=./log/master_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --partition=ivy,sbr,hwell,hpodall,spodall
#SBATCH --ntasks=1

# RES=$(sbatch -o $SLURM_ARRAY_TASK_ID.out calc.sh $1 $SLURM_ARRAY_TASK_ID) && echo ${RES##* } >> jobs.id
source ~/.bashrc

conda activate ommrew
module load cuda/9.2

date

if [ ! -e "python" ]; then
 mkdir python
fi


s=$SLURM_ARRAY_TASK_ID
# $1 is the traj template, $2 is the starting index of traj
python calc.py -t $1 -s $s --psf $2 > ./python/trj_$s.out 2> ./python/trj_$s.err

date
