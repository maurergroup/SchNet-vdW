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
from schnetpack.interfaces import SpkCalculator
from schnetpack.environment import SimpleEnvironmentProvider
from ase.calculators.calculator import Calculator, all_changes
from schnetpack import Properties
import numpy as np

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


class SpkVdwCalculatorError(Exception):
    pass #TODO RJM This doesnt look right. Errors should never pass? Is this being used?


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
            **kwargs
    ):

        #query_results = {}
        #initialise base class
        self.nmodels = nmodels
        #number of hirshfeld models
        self.nhmodels = nhmodels
        # idea is to initialize n calculators that all have a different model
        QueryCalculator.__init__(
                    self,
                    model,#[qbc_model],
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
                    *kwargs)
      

        #do additional things that are not already done by base init
        self.hirshfeld_model = hirshfeld_model
        #self.model = model
        self.models = model
        self.energy_shift = energy_shift
        self.fmax = fmax
        if self.hirshfeld_model is not None:
            for qbc_h in range(self.nhmodels):
                self.hirshfeld_model[qbc_h] = self.hirshfeld_model[qbc_h].to(device)
            self.hirsh_volrat=hirsh_volrat
            self.model_hirshfeld = hirsh_volrat

                  
            
        #adapt the fmax value
        #this is the default
        self.adaptive_fmax = self.fmax
        def get_hirsh_volrat(self): 
            if ('output' in self.parameters and
            'hirsh_volrat' not in self.parameters['output']):
                 raise NotImplementedError
            return Calculator.get_property(self, 'hirsh_volrat', self.atoms)
        def get_adaptive_fmax(self): 
            # I cannot ask for this, spk does not have any property etc.
            # but we can always return it as a factor to multiply fmax   
            return self.adaptive_fmax
    
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

        if self.calculation_required(atoms, properties):
            if self.model_hirshfeld is not None:
                self.query_results["hirsh_volrat"] = np.zeros((self.nhmodels,len(atoms)))
                model_inputs = self.atoms_converter(atoms)
                for qbc_h in range(self.nhmodels):
                    hirshfeld_model_results = self.hirshfeld_model[qbc_h](model_inputs)
                    if self.hirsh_volrat not in hirshfeld_model_results.keys():
                        raise SpkVdwCalculatorError(
                           "Your model does not support hirshfeld volume rations. Please check the model"
                           )
                    hirshfeld = hirshfeld_model_results[self.hirsh_volrat].cpu().data.numpy()
                    self.query_results["hirsh_volrat"][qbc_h]=hirshfeld.reshape(-1)
                self.results["hirsh_volrat"] = np.mean(self.query_results["hirsh_volrat"],axis=0)
                self.results["hirsh_volratmean"] = np.mean(self.query_results["hirsh_volrat"],axis=0)
                self.results["hirsh_volratvar"] = np.var(self.query_results["hirsh_volrat"],axis=0)
         
        #save the different predictions in self.query_results
        # take the mean and std later
        #print((qbc,"step")
        print("Energy mean ",self.results["energymean"])
        print("Energy var ", self.results["energyvar"])
        print("Fmax ", self.results["fmax"])
        print("Fmax mean", self.results["fmaxmean"])
        print("Forcesmax var", self.results["fmaxvar"])
        if "hirsh_volrat" in self.results:
            print("Hirsh mean", self.results["hirsh_volratmean"])
            print("Hirsh var", self.results["hirsh_volratvar"])
        if "energystd" in self.results:
            # JW we want 0.05 and can add the error of the models as "noise" ? forces are tricky as I don't know which indices due to constraints. 
            # for instance, some atoms are always fixed and there are super large forces 
            # the stds are in the range of 0.001-0.005, which fits quite well in an interpolative regime
            self.adaptive_fmax = self.adaptive_fmax + self.results["energystd"]

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
    implemented_properties = [energy, forces, stress]
    def __init__(
        self,
        model,
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
        **kwargs
    ):
        self.nmodels = nmodels
        self.Calculators = {}
        self.query_results = {}
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
                        raise SpkCalculatorError(
                            "'{}' is not a property of your model. Please "
                            "check the model "
                            "properties!".format(self.model_energy)
                        )
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
                        raise SpkCalculatorError(
                            "'{}' is not a property of your model. Please "
                            "check the model"
                            "properties!".format(self.model_forces)
                        )
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
                        raise SpkCalculatorError(
                            "Cell with 0 volume encountered for stress computation"
                        )

                    if self.model_stress not in model_results.keys():
                        raise SpkCalculatorError(
                            "'{}' is not a property of your model. Please "
                            "check the model"
                            "properties! If desired, stress tensor computation can be "
                            "activated via schnetpack.utils.activate_stress_computation "
                            "at ones own risk.".format(self.model_stress)
                        )
                    stress = model_results[self.model_stress].cpu().data.numpy()
                    self.query_results[self.stress] = stress.reshape((3, 3)) * self.stress_units

                self.query_results = self.query_results
        for prop in self.query_results:
            if prop == "forces": # in  self.query_results:
                self.results[prop] = np.mean(self.query_results[prop],axis=0)
                self.results[prop+"mean"] = np.mean(self.query_results[prop],axis=0)
                self.results[prop+"var"] = np.var(self.query_results[prop],axis=0)
                self.results["fmaxmean"] = np.max(np.abs(self.results["forcesmean"]))
                fmax = np.max(np.abs(self.results["forces"]),axis=0)
                self.results["fmax"] = np.max(fmax)
                self.results["fmaxvar"] = np.var(fmax)
            else:
                self.results[prop] = np.mean(self.query_results[prop])
                self.results[prop+"mean"] = np.mean(self.query_results[prop])
                self.results[prop+"var"] = np.var(self.query_results[prop])
            
 
