#!/bin/bash
#SBATCH --job-name=fchk-rdf
#SBATCH --time=00:10:00
#SBATCH --partition=sbr
#SBATCH --ntasks=1

# initialize job id array
echo 1 > done.rdf
