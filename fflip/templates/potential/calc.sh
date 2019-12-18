#!/bin/bash
#SBATCH --output=./log/potc_%A_%a.out
#SBATCH --error=./log/potc_%A_%a.err
#SBATCH --array=0-43
#SBATCH --time=03:00:00
#SBATCH --partition=ivy,sbr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=3G

source ~/.bashrc

conda activate ommrew-7.4
module load cuda/10.1

date

#if [ ! -e "slurmq_job_array_files" ]; then
# mkdir slurmq_job_array_files
#fi

if [ ! -e "python" ]; then
 mkdir python
fi

if [ ! -e "block_data" ]; then
 mkdir block_data
fi

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
python calc.py -t $1 -s $2 --psf $3 --percentage $4 > ./python/trj_$2_$SLURM_ARRAY_TASK_ID.out 2> ./python/trj_$2_$SLURM_ARRAY_TASK_ID.err

date
