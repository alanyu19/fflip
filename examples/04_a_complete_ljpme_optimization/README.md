# Nonbonded Term Optimization

Workflow to optimize nonbonded parameters in C36 lipid force field

<h3> Indexing of properties in this example:</h3>

0: area dppc bilayer @323.15K  
1: area dppc bilayer @333.15K  
2. overall thickness dppc bilayer @323.15K
3. overall thickness dppc bilayer @333.15K  

<h3> Do these upates:</h3>

step 1. Define your lipids in **$path_to_fflip/fflip/chm/lipids** and import it in **$path_to_fflip/fflip/chm/util.py** (from fflip.chm.lipids.$your_lipid import \*)  
step 2. In **$path_to_fflip/fflip/ljpme/scheme.py**, find the function called *lipfinder* and add your lipid with a preferred name (note this "name" will be used in targets.py)

*after doing step 1 & step 2, run **setup.py*** in **$path_to_fflip** *to install fflip (I recommend using "python setup.py develop")*   
`* some of the required packages will not be be installed automatically, please read the main README for more info` 

step 3. Modify the simulation template in the **templates** folder (in most cases, only change to sdyn.sh is needed)  
step 4. Create your template(s) for property calculation in the **templates** folder (take a look at the **area** folder!)  
step 5. Modify the potential calculation template in the **templates** folder  
step 6. Put your psf files and crd files in **psf_files** and **crd_files**, respectively  
step 7. Update the paths in **targets.py** (see comments in this file for instructions)  
step 8. Update your targets in **targets.py**

<h2> After these:</h2>

(1) Run your simulation by: ***fflip simulate [OPTIONS]*** (use "--help" to see what you can do with this)  
(2) Run your observable/potential calculations by ***fflip obsopt [OPTIONS]***  
(3) Run your reweighting by ***fflip scalc [OPTIONS]***  
(4) Run your error analysis on the sensitivity by ***fflip rcalc [OPTIONS]***  
(5) Run a linear optimization by ***fflip linearopt [OPTIONS]***
