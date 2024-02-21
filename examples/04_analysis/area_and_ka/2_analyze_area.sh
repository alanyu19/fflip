# Note: -s is the time interval between your datapoints in area.dat, which depends on your setting in dyn.py
fflip analyze-area-per-lipid -i area.dat -l 36 -t 323.15 -b 10 -s 0.004


# FFLiP: A Python Package for CHARMM/Drude Lipid FF Parameterization
# Usage: fflip analyze-area-per-lipid [OPTIONS]
# 
#   Equilibrium detection and area per lipid plot using area per lipid time
#   series data
# 
# Options:
#   -i, --input_file TEXT        Input area data file
#   -l, --nlip INTEGER           Number of lipids per leaflet
#   -t, --temperature FLOAT      Temperature of the simulation (Kelvin)
#   -b, --block_size FLOAT       Block size to use for statistical analysis,
#                                for phospholipids 10 ns is recommended (also is
#                                default)
#   -s, --interval FLOAT         Time interval (ns) between data points
#   -f, --force                  If this flag provided, force the calculation
#                                even if not sufficient equlibrium data
#   -ff, --force_fraction FLOAT  Proportion of data forced to use (ignoring the
#                                equilibrium starting point)
#   --help                       Show this message and exit.
