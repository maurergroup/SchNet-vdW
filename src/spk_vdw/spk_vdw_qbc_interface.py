"""
This module provides a an ASE calculator class [#ase1]_ for SchNetPack models derived from SpkCalculator. 
It adds the option to predict hirshfeld ratios for vdW corrections.

References
----------
.. [#ase1] Larsen, Mortensen, Blomqvist, Castelli, Christensen, Dułak, Friis,
    Groves, Hammer, Hargus: The atomic simulation environment -- a Python
    library for working with atoms.
    Journal of Physics: Condensed Matter, 9, 27. 2017.

.. [#spk1] 
"""

import os
import torch
from ase import units
import ase.io
from schnetpack.interfaces import SpkCalculator
from schnetpack.environment import SimpleEnvironmentProvider
from ase.calculators.calculator import Calculator, all_changes
from schnetpack import Properties
import numpy as np
import json
from schnetpack.data.atoms import AtomsConverter
from schnetpack.utils.spk_utils import DeprecationHelper
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.io.xyz import read_xyz, write_xyz
from ase.md import VelocityVerlet, Langevin, MDLogger
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations

from schnetpack.md.utils import MDUnits

from schnetpack.environment import SimpleEnvironmentProvider




class SpkVdwCalculator(SpkCalculator):
    """
    ASE calculator for schnetpack machine learning models with vdW.

    Args:
        ml_model (schnetpack.AtomisticModel): Trained model for
            calculations
        device (str): select to run calculations on 'cuda' or 'cpu'
        collect_triples (bool): Set to True if angular features are needed,
            for example, while using 'wascf' models
        environment_provider (callable): Provides neighbor lists
        pair_provider (callable): Provides list of neighbor pairs. Only
            required if angular descriptors are used. Default is none.
        **kwargs: Additional arguments for basic ase calculator class
    """
    energy = Properties.energy
    forces = Properties.forces
    stress = Properties.stress
    hirsh_volrat = "hirshfeld_volumes"
    implemented_properties = [energy, forces, stress, "hirsh_volrat"]
    def __init__(
            self,
            model,
            hirshfeld_model = None,
            device="cpu",
            collect_triples=False,
            environment_provider=SimpleEnvironmentProvider(),
            energy=None,
            forces=None,
            stress=None,
            hirsh_volrat=None,
            energy_units="eV",
            forces_units="eV/Angstrom",
            stress_units="eV/Angstrom/Angstrom/Angstrom",
            energy_shift=None,
            fmax = 0.05,
            nmodels = int(1),
            nhmodels = int(1),
            extrapolate = False,
            adaptive_fmax = None,
            qbc = None,
            logfile = 'opt.log',
            dftd4 = None,
            functional = "pbe",
            **kwargs
    ):

        #query_results = {}
        #initialise base class
        self.qbc = qbc
        self.dftd4 = dftd4
        self.nmodels = nmodels
        #number of hirshfeld models
        self.nhmodels = nhmodels
        self.functional = functional
        # idea is to initialize n calculators that all have a different model
        QueryCalculator.__init__(
                    self,
                    model,#[qbc_model],
                    hirshfeld_model,
                    device,
                    collect_triples,
                    environment_provider,#[qbc],
                    energy,
                    forces,
                    stress,
                    energy_units=energy_units,
                    forces_units=forces_units,
                    stress_units=stress_units,
                    energy_shift=energy_shift,
                    nmodels=nmodels,
                    nhmodels=nhmodels,
                    extrapolate=extrapolate,
                    adaptive_fmax = adaptive_fmax,
                    qbc = qbc,
                    logfile = logfile,
                    **kwargs)
      
        #do additional things that are not already done by base init
        self.hirshfeld_model = hirshfeld_model
        #self.model = model
        self.models = model
        self.energy_shift = energy_shift
        self.fmax = fmax
        self.logfile = logfile
        if self.hirshfeld_model is not None:
            for qbc_h in range(self.nhmodels):
                self.hirshfeld_model[qbc_h] = self.hirshfeld_model[qbc_h].to(device)
            self.hirsh_volrat=hirsh_volrat
            self.model_hirshfeld = hirsh_volrat

                  
        # define if system is known to the model    
        # false if structure is included in the training set; false is default
        self.bh = kwargs["mode"]
        self.extrapolate = extrapolate
        self.current_fmax = np.inf
        self.previous_fmax = self.current_fmax
        self.reached_fmax = False
        self.last_evar = np.inf 
        self.vdwener = np.inf
        self.nsteps = int(0)
        self.end = False
        self.adaptive_fmax = adaptive_fmax
        def get_hirsh_volrat(self): 
            if ('output' in self.parameters and
            'hirsh_volrat' not in self.parameters['output']):
                 raise NotImplementedError
            return Calculator.get_property(self, 'hirsh_volrat', self.atoms)

        """def get_end_of_optimization(self, current_fmax):  
            init_=(0,np.inf)
            print(init_)
            self.current_fmax = current_fmax
            np.savez("qbc_evar.npz",init_)
            print(self.current_fmax)
            return Calculator.get_propy(self, 'end', self.current_fmax)"""
                                           
    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        """                                
        Args:                              
            atoms (ase.Atoms): ASE atoms object.
            properties (list of str): do not use this, no functionality
            system_changes (list of str): List of changes for ASE.
        """
 
        QueryCalculator.calculate(
            self,
            atoms,
            properties, 
            system_changes
        )
        # save the previous step
        if self.dftd4 == True:
            # read latest geometry step
            try:
                # for the first step, tmp.xyz is already provided
                latest_geom = ase.io.read(self.logfile[:-3]+"traj",-1)
                ase.io.write("tmp.xyz",latest_geom)
            except:
                pass
            # write out latest geometry step
            # conduct dftd4 calculation
            os.system("dftd4 --func %s --json --grad --noedisp tmp.xyz" %self.functional)

            # read in json file that is written by dftd4
            f = open('dftd4.json')
            dftd4data = json.load(f)
            print(dftd4data)
            print(self.results.keys())
            self.results["energy"] -= np.array(dftd4data["energy"])
            # + because gradient
            self.results["forces"] += np.array(dftd4data["gradient"]).reshape(int(len(dftd4data["gradient"])/3),3)
        if self.qbc != None: 
          try:
            file = open(self.logfile,"r").readlines()    
            if len(file) >= 3:
                if file[len(file)-1].startswith("BasinHopping"):
                    self.current_fmax = np.inf 
                    self.last_evar = np.inf
                    self.last_fmax = np.inf
                    self.nsteps = 0
                    self.nsteps_fmax = 0
                    self.bh = True
                    bhstep = file[len(file)-1].split()[2]
                    #os.system("mv Stopped Stopped_%s"%bhstep[:len(bhstep)])
                    self.previous_fmax = np.inf
                    #next bh step starts, set fmax reached to False
                    self.reached_fmax=False 
                    self.end = False
                    self.previous_ener = self.vdwener
                elif file[len(file)-1].startswith("BFGS") or  file[len(file)-1].startswith("LBFGS") :
                    if self.current_fmax:
                        self.previous_fmax = self.current_fmax #float(file[len(file)-1].split()[
                        self.previous_ener = self.vdwener
                    else:
                        self.previous_fmax = np.inf
                    self.current_fmax = float(file[len(file)-1].split()[4])
                    self.vdwener = float(file[len(file)-1].split()[3])
                    if self.previous_ener == np.inf:
                        #self.previous_ener = float(file[2].split()[3])
                        self.previous_ener = self.vdwener
                    if self.reached_fmax == True:
                        self.last_evar = float(file[len(file)-2].split()[2])
                        self.nsteps = int(file[len(file)-2].split()[1])
                else:
                    pass
            else:
                self.current_fmax = 1.9
                self.previous_fmax = 1.9
                self.last_evar = np.inf
                self.nsteps = 0
                self.nsteps_fmax = 0
                self.previous_ener = np.inf
                self.vdwener = np.inf
          except IOError:
            self.current_fmax = np.inf 
            self.last_evar = np.inf
            self.nsteps = 0
            self.nsteps_fmax = 0
            self.previous_ener = np.inf 
            #print("Please save the output in a file named %s (python ..py >> %s) for using the QBC models and early stopping."%(self.logfile,self.logfile))
          if self.extrapolate == False:
            if self.current_fmax < self.adaptive_fmax or  self.reached_fmax == True or (self.current_fmax > 2*self.previous_fmax and self.bh==False) or (self.current_fmax > 3 and self.bh==False) or abs(self.vdwener-self.previous_ener)>1:
                self.reached_fmax = True
                if self.results['energyvar'] > self.last_evar:
                    self.nsteps+=1
                    
                    print("Evar: %i %f" %(self.nsteps,self.results['energyvar']))
                else:
                    self.nsteps = 0
                    print("Evar: %i %f " %(self.nsteps,self.results['energyvar']))
                if self.nsteps >= int(4)  or (abs(self.vdwener-self.previous_ener)>1 and self.bh==False) or (self.current_fmax>3 and self.bh==False) or (self.results["energyvar"]>1 and self.bh==False):
                    self.end = True
              
                if self.current_fmax > self.previous_fmax:
                    self.nsteps_fmax += 1
                else:
                    self.nsteps_fmax = 0
                if self.nsteps_fmax >= int(4):
                    self.end = True
            else:
                if self.current_fmax > self.previous_fmax:
                    self.nsteps_fmax += 1
                else:
                    self.nsteps_fmax = 0
                if self.nsteps_fmax >= int(4):
                    self.end = True
          else:
            if self.reached_fmax == True:
                self.end = True
                self.reached_fmax = True
                if self.results['energyvar']> self.last_evar:
                    self.end = False
                    print("Evar: %i %f" %(self.nsteps,self.results['energyvar']))
                else: 
                    print("Evar: %i %f" %(self.nsteps,self.results['energyvar']))
                    if self.current_fmax > self.previous_fmax:
                        self.nsteps_fmax += 1
                    else:
                        self.nsteps_fmax = 0
                    if self.nsteps_fmax >= int(3):
                        self.end = True
                    else: 
                        self.end = False
 
          if self.end == True:
            if self.results['energyvar']>1: 
                os.system("echo 'Terminated optimization early with energy variance %f.' >> Stopped_untrustworthy"%self.results['energyvar'])
                #os.system("mv Stopped_untrustworthy Stopped_untrusstworthy_%s"%bhstep[:len(bhstep)])
            else: 
                pass
                #os.system("echo 'Terminated optimization early but still likely to be trustworthy.' >> Stopped")
          if self.end == True and self.bh == False:
                self.results["forces"]=self.results["forces"]*0.0 #exit()
                self.results["forces"]=np.ma.masked_equal(self.results["forces"],0.0)
          """print("Energy mean ",self.results["energymean"])
          print("Energy var ", self.results["energyvar"])
          print("Fmax ", self.results["fmax"])
          print("Fmax mean", self.results["fmaxmean"])
          print("Forcesmax var", self.results["fmaxvar"])
          if "hirsh_volrat" in self.results:
              print("Hirsh var", np.mean(self.results["hirsh_volratvar"]))
          if "energystd" in self.results:
              # JW we want 0.05 and can add the error of the models as "noise" ? forces are tricky as I don't know which indices due to constraints. 
              # for instance, some atoms are always fixed and there are super large forces 
              # the stds are in the range of 0.001-0.005, which fits quite well in an interpolative regime
              self.adaptive_fmax = self.adaptive_fmax + self.results["energystd"]"""

class QueryCalculator(Calculator):
    """
    ASE calculator for schnetpack machine learning models.

    Args:
        ml_model (schnetpack.AtomisticModel): Trained model for
            calculations
        device (str): select to run calculations on 'cuda' or 'cpu'
        collect_triples (bool): Set to True if angular features are needed,
            for example, while using 'wascf' models
        environment_provider (callable): Provides neighbor lists
        pair_provider (callable): Provides list of neighbor pairs. Only
            required if angular descriptors are used. Default is none.
        **kwargs: Additional arguments for basic ase calculator class
    """

    energy = Properties.energy
    forces = Properties.forces
    stress = Properties.stress
    implemented_properties = [energy, forces, stress, "hirsh_volrat"]
    def __init__(
        self,
        model,
        hirshfeld_model,
        device="cpu",
        collect_triples=False,
        environment_provider=SimpleEnvironmentProvider(),
        energy=None,
        forces=None,
        stress=None,
        energy_units="eV",
        forces_units="eV/Angstrom",
        stress_units="eV/Angstrom/Angstrom/Angstrom",
        nmodels = int(1),
        nhmodels = int(1),
        **kwargs
    ):
        self.nmodels = nmodels
        self.nhmodels = nhmodels
        self.Calculators = {}
        self.query_results = {}
        self.results = {}
        self.model_hirshfeld = None
        for qbc in range(self.nmodels):
            Calculator.__init__(self, **kwargs)
            self.Calculators[qbc] = Calculator
            self.models = model
            #self.model.to(device)
            self.device = device
            self.atoms_converter = AtomsConverter(
                environment_provider=environment_provider[qbc],
                collect_triples=collect_triples,
                device=device,
            )
            self.model_energy = energy
            self.model_forces = forces
            self.results = {}
            self.model_stress = stress
            # Convert to ASE internal units (energy=eV, length=A)
            self.energy_units = MDUnits.unit2unit(energy_units, "eV")
            self.forces_units = MDUnits.unit2unit(forces_units, "eV/Angstrom")
            self.stress_units = MDUnits.unit2unit(stress_units, "eV/A/A/A")
        for qbc_h in range(self.nhmodels):
            Calculator.__init__(self, **kwargs)
            self.Calculators[self.nmodels+qbc_h] = Calculator
            self.hmodels = hirshfeld_model
            #self.model.to(device)
            self.device = device
            self.atoms_converter = AtomsConverter(
                environment_provider=environment_provider[qbc],
                collect_triples=collect_triples,
                device=device,
            )
            self.model_hirshfeld = self.hirsh_volrat
            # Convert to ASE internal units (energy=eV, length=A)
 
    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        """
        Args:
            atoms (ase.Atoms): ASE atoms object.
            properties (list of str): do not use this, no functionality
            system_changes (list of str): List of changes for ASE.
        """
        # First call original calculator to set atoms attribute
        # (see https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/calculator.html#Calculator)
        if self.calculation_required(atoms, properties):
            for qbc in range(self.nmodels):
                self.Calculators[qbc].calculate(self, atoms)
                # Convert to schnetpack input format
                model_inputs = self.atoms_converter(atoms)
                # Call model
                self.model = self.models[qbc].to(self.device)
                model_results = self.model(model_inputs)
                # Convert outputs to calculator format
                if self.model_energy is not None:
                    if self.model_energy not in model_results.keys():
                        print("'{}' is not a property of your model. Please check the model properties!".format(self.model_energy))
                    energy = model_results[self.model_energy].cpu().data.numpy()
                    if self.energy in self.query_results:
                        self.query_results[self.energy][qbc]=(energy.item() * self.energy_units)
                    else:
                        self.query_results[self.energy] = np.zeros((self.nmodels,1))
                        self.query_results[self.energy][qbc]  = (
                            energy.item() * self.energy_units
                        )  # ase calculator should return scalar energy

                if self.model_forces is not None:
                    if self.model_forces not in model_results.keys():
                        print("'{}' is not a property of your model. Please check the model properties!".format(self.model_forces))
                    forces = model_results[self.model_forces].cpu().data.numpy()
                    if self.forces in self.query_results:
                        self.query_results[self.forces][qbc] = (
                        forces.reshape((len(atoms), 3)) * self.forces_units
                        )
                    else:
                        self.query_results[self.forces] = np.zeros((self.nmodels,len(atoms),3))
                        self.query_results[self.forces][qbc]= (forces.reshape((len(atoms),3)) * self.forces_units )

                if self.model_stress is not None:
                    if atoms.cell.volume <= 0.0:
                        print("Cell with 0 volume encountered for stress computation")

                    if self.model_stress not in model_results.keys():
                        print("'{}' is not a property of your model. Please check the model properties! Set to 0 for now".format(self.model_stress))
                        stress = model_results[self.model_stress] = np.zeros((3,3))
                    else:
                        stress = model_results[self.model_stress].cpu().data.numpy()
                    self.query_results[self.stress] = stress.reshape((3, 3)) * self.stress_units
                
            for qbc_h in range(self.nhmodels):
                if self.model_hirshfeld is not None:
                    self.Calculators[self.nmodels + qbc_h].calculate(self, atoms)
                    # Convert to schnetpack input format
                    model_inputs = self.atoms_converter(atoms)
                    # Call model
                    self.model = self.hmodels[qbc_h].to(self.device)
                    model_results = self.model(model_inputs)
                    if "hirsh_volrat" not in self.query_results:
                        self.query_results["hirsh_volrat"] = np.zeros((self.nhmodels,len(atoms)))
                    if self.hirsh_volrat not in model_results.keys():
                        raise SpkVdwCalculatorError(
                           "Your model does not support hirshfeld volume rations. Please check the model"
                           )
                    hirshfeld = model_results[self.hirsh_volrat].cpu().data.numpy()
                    self.query_results["hirsh_volrat"][qbc_h]=hirshfeld.reshape(-1)

            self.query_results = self.query_results
        for prop in self.query_results:
            if prop == "forces": # in  self.query_results:
                # test to take only one model nr. 1
                self.results[prop] = np.mean(self.query_results[prop],axis=0) #= self.query_results[prop][0] #
                self.results[prop+"mean"] = np.mean(self.query_results[prop],axis=0)
                self.results[prop+"var"] = np.var(self.query_results[prop],axis=0)
                self.results["fmaxmean"] = np.max(np.abs(self.results["forcesmean"]))
                fmax = np.max(np.abs(self.results["forces"]),axis=0)
                self.results["fmax"] = np.max(fmax)
                self.results["fmaxvar"] = np.var(fmax)
            elif prop=="hirsh_volrat":
                self.results["hirsh_volrat"] = np.mean(self.query_results["hirsh_volrat"],axis=0)
                self.results["hirsh_volratmean"] = np.mean(self.query_results["hirsh_volrat"],axis=0)
                self.results["hirsh_volratvar"] = np.var(self.query_results["hirsh_volrat"],axis=0)
            else:
                self.results[prop] = np.mean(self.query_results[prop])
                self.results[prop+"mean"] = np.mean(self.query_results[prop])
                self.results[prop+"var"] = np.var(self.query_results[prop])
              

