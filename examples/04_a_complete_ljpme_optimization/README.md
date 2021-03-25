# Nonbonded Term Optimization

Workflow to optimize nonbonded parameters in the C36 lipid force field

<h3> Indexing of properties in this example:</h3>

0: area dppc bilayer @323.15K  
1: area dppc bilayer @333.15K  
2. overall thickness dppc bilayer @323.15K  
3. overall thickness dppc bilayer @333.15K  

<h3> Do these updates:</h3>

step 1. Define your lipid(s) in **def_lipids.py** and make sure to add your lipid(s) to "match_lipid" at the end of the file  
step 2. Modify the simulation template in the **templates** folder (in most cases, only change to sdyn.sh is needed)  
step 3. Create your template(s) for property calculation in the **templates** folder (take a look at the **area** folder!)  
step 4. Modify the potential calculation template in the **templates** folder  
step 5. Copy your **def_lipids.py** to your simulation template folder and your potential calculation template folder  
step 6. Put your psf files and crd files in **psf_files** and **crd_files**, respectively  
step 7. Update the paths in **targets.py** (see comments in this file for instructions)  
step 8. Update your targets in **targets.py**

<h2> After these:</h2>

(1) Run your simulation by: ***fflip simulate [OPTIONS]*** (use "--help" to see what you can do with this)  
(2) Run your observable/potential calculations by ***fflip obsopt [OPTIONS]***  
&nbsp;&nbsp;&nbsp;*in this example, the following commands were used:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fflip obspot -i 0,1,2,3 -tl simulation -f 11 -l 40 -c observable -p 0.1 --iteration 0 (4 observables)   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fflip obspot -i 0,1 -tl simulation -f 11 -l 40 -c potential -p 0.1 --iteration 0 (only 2 systems)*  

(3) Run your reweighting by ***fflip scalc [OPTIONS]***  
&nbsp;&nbsp;&nbsp;*in this example, the following command was used:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fflip scalc -i all -f 11 -l 40 -p 0.1 --iteration 0 -P hwell*  

(4) Run your error analysis on the sensitivity by ***fflip rcalc [OPTIONS]***  
&nbsp;&nbsp;&nbsp;*in this example, the following command was used:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fflip rcalc -i all -f 11 -l 40 -p 0.1 --iteration 0 -P hwell*  

(5) Run a linear optimization by ***fflip linearopt [OPTIONS]***
&nbsp;&nbsp;&nbsp;*in this example, the following command was used:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fflip linearopt -i 0 -p 0.1 -s 0.06 -e 0.06 -c 0.03 --ssr*  

<h4> Now you can find FFLiP's solution for nonbonded parameters in the solutions folder!<\h4>  

(6) Repeat 1-5 for more iterations (some OPTIONs should be different from the zeroth iteration, use "--help" to find out!)
