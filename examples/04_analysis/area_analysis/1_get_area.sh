# Note: the dyn{}.dcd format is the standard output from fflip, which uses rickflow/openmm as the backend.
# Curly brace is the standard formatter for python.
# UPDATE YOUR PATHs ACCORDINGLY!
fflip area-per-lipid -p ../sim/dppc_ua.psf -l 36 -t ../sim/trj/dyn{}.dcd -f 1


# Output of fflip area-per-lipid --help:

# FFLiP: A Python Package for CHARMM/Drude Lipid FF Parameterization
# Usage: fflip area-per-lipid [OPTIONS]
# 
#   Area per Lipid
# 
# Options:
#   -p, --psf_file TEXT        PSF file used to run the simulation
#   -l, --nlip INTEGER         Number of lipids. Default is 36.
#   -t, --traj_template TEXT   Trajetory file template, use curly brace for
#                              index, default is "trj/dyn{}.dcd"
#   -f, --first_seqno INTEGER  first seqno to calculate the area
#   -e, --last_seqno INTEGER   last seqno to calculate the area (optional)
#   --help                     Show this message and exit.
