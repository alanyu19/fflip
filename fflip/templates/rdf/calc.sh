#!/bin/bash
#SBATCH --job-name=rdf-calc
#SBATCH --output=./log/rdfc_%A_%a.out
#SBATCH --error=./log/rdfc_%A_%a.err
#SBATCH --array=1-8
#SBATCH --time=06:00:00
#SBATCH --partition=ivy,sbr,hwell,spodall
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G

source ~/.bashrc

conda activate ommrew-7.4
module load cuda/10.1

date


if [ ! -e "python" ]; then
 mkdir python
fi

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
# $SLURM_ARRAY_TASK_ID is used for selecting the rdf pair(s)
# which are defined in calc.py
python calc.py -t $1 -s $2 --psf $3 > ./python/trj_$2_$SLURM_ARRAY_TASK_ID.out 2> ./python/trj_$2_$SLURM_ARRAY_TASK_ID.err

date
