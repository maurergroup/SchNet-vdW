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
            **kwargs
    ):

        #query_results = {}
        #initialise base class
        self.query_results = {}
        self.nmodels = nmodels
        #not sure if the loop is needed here
        for qbc in range(nmodels):
            SpkCalculator.__init__(
                    self,
                    model[qbc],#[qbc_model],
                    device,
                    collect_triples,
                    environment_provider[qbc],
                    energy,
                    forces,
                    stress,
                    energy_units=energy_units,
                    forces_units=forces_units,
                    stress_units=stress_units,
                    energy_shift=energy_shift,
               *kwargs)
      

            #do additional things that are not already done by base init
            self.hirshfeld_model = hirshfeld_model[qbc]
            self.model = model[qbc]
            self.energy_shift = energy_shift
            self.fmax = fmax
            if self.hirshfeld_model is not None:
                for qbc in range(nmodels):
                    self.hirshfeld_model = self.hirshfeld_model.to(device)
            else:
                self.hirsh_volrat=hirsh_volrat
                self.model_hirshfeld = hirsh_volrat

                  
            
        #adapt the fmax value
        #this is the default
        self.qbc = qbc
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
        for qbc in range(self.nmodels):
            self.qbc=qbc
            SpkCalculator.calculate(
                self,
                atoms,
                properties, 
                system_changes
            )

            if self.calculation_required(atoms, properties):
                if self.model_hirshfeld is not None:
                    model_inputs = self.atoms_converter(atoms)
                    hirshfeld_model_results = self.hirshfeld_model(model_inputs)
                    if self.hirsh_volrat not in hirshfeld_model_results.keys():
                        raise SpkVdwCalculatorError(
                           "Your model does not support hirshfeld volume rations. Please check the model"
                           )
                    hirshfeld = hirshfeld_model_results[self.hirsh_volrat].cpu().data.numpy()
                    self.results["hirsh_volrat"]=hirshfeld.reshape(-1)
            #save the different predictions in self.query_results
            # take the mean and std later
            #print(self.qbc,"step")
            if "energy" not in self.query_results and "energy" in self.results:
                 if self.energy_shift is not None:
                     self.results["energy"] += self.energy_shift
                 self.query_results["energy"]=np.zeros((self.nmodels,1))
            self.query_results["energy"][self.qbc]=self.results["energy"]

            if "forces" not in self.query_results and "forces" in self.results:
                self.query_results["forces"]=np.zeros((self.nmodels,len(self.results["forces"]),3))
            self.query_results["forces"][self.qbc]=self.results["forces"]
            if "hirsh_volrat" not in self.query_results and self.hirsh_volrat in self.results:#self.hirshfeld_model_results:
                self.query_results["hirsh_volrat"] = np.zeros((self.nmodels,len(self.hirshfeld_model_results[self.hirsh_volrat])))
            if "hirsh_volrat" in self.results:#self.hirshfeld_model_results:
                self.query_results["hirsh_volrat"][qbc] = self.results[hirsh_volrat]
            for prop in self.query_results:
                self.results[prop] = np.mean(self.query_results[prop],axis=0)
                self.results[prop+"std"] = np.std(self.query_results[prop])
            if "energystd" in self.results:
                # JW we want 0.05 and can add the error of the models as "noise" ? forces are tricky as I don't know which indices due to constraints. 
                # for instance, some atoms are always fixed and there are super large forces 
                # the stds are in the range of 0.001-0.005, which fits quite well in an interpolative regime
                self.adaptive_fmax = self.adaptive_fmax + self.results["energystd"]
