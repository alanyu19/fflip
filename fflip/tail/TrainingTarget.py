#!/usr/bin/env python

from __future__ import division, print_function


import numpy as np
import MDAnalysis as mda
import pandas as pd

import os
import time
import sys


def bzavg(obs,boltz):
    # Get the Boltzmann average of an observable.
    if obs.ndim == 2:          # ndim is the dimension of data set, a method of numpy
        if obs.shape[0] == len(boltz) and obs.shape[1] == len(boltz):   # number of rows and columns """
            raise Exception('Error - both dimensions have length equal to number of snapshots, now confused!')
        elif obs.shape[0] == len(boltz):
            return np.sum(obs*boltz.reshape(-1,1),axis=0)/np.sum(boltz)  # now the average is got
        elif obs.shape[1] == len(boltz):
            return np.sum(obs*boltz,axis=1)/np.sum(boltz)
        else:
            raise Exception('The dimensions are wrong!')  # when none of them has the correct length
    elif obs.ndim == 1:
        return np.dot(obs,boltz)/sum(boltz)
    else:
        raise Exception('The number of dimensions can only be 1 or 2!')

def calc_kappa(b=None, **kwargs):
    #global kT
    #global nframes
    if b is None: b = np.ones(kwargs["nframes"], dtype=float)
    if 'v_' in kwargs:
        v_ = kwargs['v_']
        # return bar_unit / kT * (bzavg(v_**2,b)-bzavg(v_,b)**2)/bzavg(v_,b)
        return (bzavg(v_**2,b)-bzavg(v_,b)**2)/bzavg(v_,b)*1.0e-30/ kwargs["kT"]
    
class TrainingTarget(object):
    
    def __init__(self, **kwargs):
        assert "name" in kwargs        
        self.name = kwargs["name"]
        assert "temp" in kwargs
        self.temp = kwargs["temp"]
        assert "resname" in kwargs
        self.resname = kwargs["resname"] 
        assert "crd" in kwargs
        self.crd = os.path.abspath(kwargs["crd"])
        assert "guess_box" in kwargs
        self.guess_box = kwargs["guess_box"]
        assert "hovnb" in kwargs
        self.hovnb = kwargs["hovnb"]
        assert "fftx" in kwargs
        self.fftx = kwargs["fftx"] 
        assert "root_dir" in kwargs
        self.dir = os.path.abspath(kwargs["root_dir"] + self.name + "_" + str(self.temp))
        assert "psf" in kwargs
        self.psf = os.path.abspath(kwargs["psf"])
        assert "number_molecules" in kwargs
        self.number = kwargs["number_molecules"]
        assert "molmass" in kwargs
        self.molmass = kwargs["molmass"]
        assert "dcd" in kwargs
        self.dcd = kwargs["dcd"]
        assert "first_dcd" in kwargs
        self.first_dcd = kwargs["first_dcd"]
        assert "last_dcd" in kwargs
        self.last_dcd = kwargs["last_dcd"]
        assert "replica" in kwargs
        self.replica = kwargs["replica"]
        assert "prop" in kwargs
        self.prop = kwargs["prop"]
        assert "exp_dic" in kwargs
        self.exp_dic = kwargs["exp_dic"]
        assert "addcalc" in kwargs
        self.addcalc = kwargs["addcalc"]
        assert "weights" in kwargs
        self.weights = kwargs["weights"]
        self.dic = {}  # Ready for filling up
        self.frames = None
        
    def __repr__(self):
        return "complete later!" 
    
    def CleanUpPreviousIteration(self):
        os.chdir(self.dir)
        os.system("echo {} > last.seqno".format(self.dcd))
        os.system("rm -rf kap.dat rho.dat *.out test.file start.time next.seqno done.* slurm-* slurm pyfile output")
        #os.system("rm -f ./output/*")
        #os.system("rm -f ./output/avg_*")
        
        if 'hov' in self.prop:
            os.system("rm -rf vap.dat gas") 
            """ done.gas ? """
        
        if 'd3d' in self.prop:
            os.system("rm -rf d3d.dat uf_trj d3d")
        
        self.dic = {} 
            
    def Simulate(self):
        os.chdir(self.dir)
        assert os.path.isfile(self.dir + "/minimize.csh")
        assert os.path.isfile(self.dir + "/mdyn.sh")
        assert os.path.isfile(self.dir + "/last.seqno")
        # pass the psf in this version since the openmm
        # rundyn.py will need that.
        # should ADD a step in minimize.inp to write out last crd
        # COPY the coordinates from minimization in mdyn.sh or rundyn.py
        os.system("sbatch minimize.csh {} {} {} {} {}".format(
                  self.temp, self.psf, self.crd, self.guess_box, self.fftx))
                
    #def CreateAnts(self):
    #    print("Creating Ants for {}".format(self.name + '_' + str(self.temp)))
    #    os.chdir(self.dir)
    #    with open("ants.trp", "w") as file:
    #        for i in range(self.number/4):
    #             file.write(str(i+1)+' ')
    
    def GetDensitykappa(self, counter):
        # Get Kappa even if no experimental data available
        os.chdir(self.dir)
        while not os.path.isfile("done.dyn"):
            time.sleep(300)
        os.chdir(self.dir)
        kT = 1.380648e-23 * float(self.temp)
        ommout = np.loadtxt("./output/dyn{}.out".format(self.last_dcd), skiprows=1, delimiter=",")
        V = ommout[:,4] * 1000
        nframes = np.shape(V)[0]
        if self.frames != None:
            assert nframes == self.frames
        else:
            self.frames = nframes
        print("Iteration #{}, Getting Density for {}".format(
              counter, self.name + '_' + str(self.temp)))
        V_avg_inverse = np.mean(1/V)    
        density = int(self.number) * float(self.molmass) * V_avg_inverse * 1e27/6.023e23
        #self.dic["rho"] = [density]
        # store the information for quick check
        assert density is not None
        os.system("echo {} > rho.dat".format(density)) 
        #assert self.dic["rho"] is not None
        print("Iteration #{}, Density for {} Ready".format(
              counter, self.name + '_' + str(self.temp)))
        with open("rho.dat") as file:
            self.dic["rho"] = [float(file.read())]
        print("Iteration #{}, Getting Kappa for {}".format(
              counter, self.name + '_' + str(self.temp)))
        kap = calc_kappa(None,**{'v_':V, 'kT':kT, 'nframes':nframes})
        os.system("echo {} > kap.dat".format(kap))
        assert kap is not None
        print("Iteration #{}, Kappa for {} Ready".format(
              counter, self.name + '_' + str(self.temp)))
        with open("kap.dat") as file:
            self.dic["kappa"] = [float(file.read())]
        sys.stdout.flush()
        
    def Unfold_and_Gas(self, counter):
        os.chdir(self.dir)
        while not os.path.isfile("done.dyn"):
            time.sleep(300)
        # Diffusion in 3D
        if 'd3d' in self.addcalc:
            #os.chdir(driver_path + system.name + "_" + str(system.temp))
            kitty = mda.Universe(self.psf, 
                                 self.dir + "/output/dyn{}.dcd".format(self.last_dcd)) 
            nframes = len(kitty.trajectory)
            if self.frames != None:
                assert nframes == self.frames
            else:
                self.frames = nframes            
            L = np.zeros(self.frames)
            for j in range(self.frames):
                L[j] = kitty.trajectory[j].triclinic_dimensions[0,0]
            L_avg = np.mean(L)
            self.mean_box_length = L_avg
            os.system("sbatch do-ReducUnfold.csh {} {} {} {} {} {}".format(
                      self.mean_box_length, self.first_dcd, self.last_dcd, self.resname,
                      self.number, self.crd))
            print("Iteration #{}, Unfolding Trajectory for {}".format(
                  counter, self.name + '_' + str(self.temp)))
            sys.stdout.flush()
        # Heat of Vaporization
        if 'hov' in self.addcalc:
            print("Iteration #{}, Creating Ants in {}".format(
                  counter, self.name + '_' + str(self.temp)))
            sys.stdout.flush()
            with open("ants.trp", "w") as file:
                for i in range(int(self.number)):
                     file.write(str(i+1)+' ')
            #os.chdir(driver_path + system.name + "_" + str(system.temp))
            os.system("sbatch gasmaster.csh {} {} {} {} {} {} {} {}".format(
                      self.last_dcd, self.number, self.temp, self.resname, 
                      self.crd, self.guess_box, self.hovnb, self.frames))
            print("Iteration #{}, Getting Gas Phase Potentials for {}".format(
                  counter, self.name + '_' + str(self.temp)))
            sys.stdout.flush()
            
    def GetHeatofVaporization(self, counter):
        if 'hov' in self.addcalc:
            os.chdir(self.dir)
            while not os.path.isfile("done.gas"):
                time.sleep(10)                
            ommout = np.loadtxt("./output/dyn{}.out".format(self.last_dcd), skiprows=1, delimiter=",")
            condensed_e = np.mean(ommout[:,2] * 0.23900573614 / self.number)
            gas_e = np.loadtxt("./gas/avg_energy.dat")          
            # R = 1.9872036(11)**10^-3 kcal/(K*mol)
            print("Iteration #{}, Gas Phase for {} Done, Calculating Heat of Vaporization".format(
                  counter, self.name + '_' + str(self.temp)))
            result = gas_e - condensed_e + self.temp * 0.0019872
            os.system("echo {} > vap.dat".format(result))
            os.system("echo 1 > done.vap")
            sys.stdout.flush()
    
    def GetDiffusionConstant(self, counter):
        # Diffusion in 3D (continued)
        if 'd3d' in self.addcalc:         
            os.chdir(self.dir)
            while not os.path.isfile("done.unf"):
                time.sleep(10)
            assert self.mean_box_length != None
            os.system("sbatch do-bcseMSD.csh {} {} {} {} {} {} {} {} {}".format(
                    self.mean_box_length, self.first_dcd, self.last_dcd,
                    self.number, self.temp, self.exp_dic["visc"],
                    self.resname, self.crd, self.frames)
                )
            print("Iteration #{}, Getting Diffusion Constant for {}".format(
                  counter, self.name + '_' + str(self.temp)))
            sys.stdout.flush()
            
    def AssertDiffusionIsDone(self, counter):
        if 'd3d' in self.addcalc:
            os.chdir(self.dir)
            while not os.path.isfile("done.msd"):
                time.sleep(10)
            with open("d3d.dat") as file:
                self.dic["d3d"] = [float(file.read())]
                assert self.dic["d3d"] is not None
                print("Iteration #{}, Diffusion Constant for {} Ready".format(
                      counter, self.name + '_' + str(self.temp)))
            sys.stdout.flush()
    
    def AsserthovIsDone(self, counter):
        if 'hov' in self.addcalc:
            os.chdir(self.dir)
            while not os.path.isfile("done.vap"):
                time.sleep(5)
            with open("vap.dat") as file:
                self.dic["hov"] = [float(file.read())]
                assert self.dic["hov"] is not None
                print("Iteration #{}, Heat of Vaporization for {} Ready".format(
                      counter, self.name + '_' + str(self.temp)))
            sys.stdout.flush()
     
    def GetSSR(self, counter):
        print("Iteration #{}, Getting SSR for {}".format(counter, self.name + '_' + str(self.temp)))
        ssr = 0
        for j, prop in enumerate(self.prop):
            ssr += self.weights[j] * ((self.dic[prop][0] - self.exp_dic[prop])/self.exp_dic[prop])**2
        self.dic["ssr"] = [ssr]
        print("Iteration #{}, SSR for {} Ready".format(counter, self.name + '_' + str(self.temp)))
        return float(self.dic["ssr"][0])
        sys.stdout.flush()
        
    def CreateTableWithExp(self):
        print("Creating New Table for {}".format(self.name + '_' + str(self.temp)))
        var_list = []
        var_list = ["CH1E_sigma", "CH1E_epsilon", "rho", "kappa"]
        for prop in self.addcalc:
            var_list.append(prop)
        var_list.append('ssr')
        var_list.append('ssr_sum')
        csv_dict = {}
        for var in var_list:
            csv_dict[var] = [self.exp_dic[var]] if var in self.exp_dic else [None]
            df = pd.DataFrame(csv_dict, columns=var_list)
            df.to_csv(self.name + '-' + str(self.temp) + ".csv", index=False)
        sys.stdout.flush()
        
    def WriteInfoToTable(self, path, counter):
        print("Iteration #{}, Writing Table for {}".format(counter, self.name + '_' + str(self.temp)))
        os.chdir(path)
        df = pd.read_csv(self.name + '-' + str(self.temp) + ".csv")
        res_df = pd.DataFrame(self.dic)
        df = pd.concat([df,res_df], sort=False)
        os.system("cp {} {}".format(self.name + '-' + str(self.temp) + ".csv",
                  self.name + '-' + str(self.temp) + "-bu.csv"))
        df.to_csv(self.name + '-' + str(self.temp) + ".csv", index=False)
        os.system("rm -f {}".format(self.name + '-' + str(self.temp)+"-bu.csv"))
        print("Iteration #{}, Table for {} Updated".format(counter, self.name + '_' + str(self.temp)))
        sys.stdout.flush()
