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

from ase import units
from schnetpack.interfaces import SpkCalculator
from schnetpack.environment import SimpleEnvironmentProvider
from ase.calculators.calculator import Calculator, all_changes
from schnetpack import Properties
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
        **kwargs
    ):

        #initialise base class
        SpkCalculator.__init__(
            self,
            model,
            device,
            collect_triples,
            environment_provider,
            energy,
            forces,
            stress,
            energy_units=energy_units,
            forces_units=forces_units,
            stress_units=stress_units,
            *kwargs)
        
        #do additional things that are not already done by base init
        self.hirshfeld_model = hirshfeld_model
        if self.hirshfeld_model is not None:
            self.hirshfeld_model.to(device)
        else:
            self.hirsh_volrat=hirsh_volrat

        self.model_hirshfeld = hirsh_volrat


    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        """
        Args:
            atoms (ase.Atoms): ASE atoms object.
            properties (list of str): do not use this, no functionality
            system_changes (list of str): List of changes for ASE.
        """

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
                self.results["hirsh_volrat"] = hirshfeld.reshape(-1)            

    def get_hirsh_volrat(self): 
        if ('output' in self.parameters and
           'hirsh_volrat' not in self.parameters['output']):
                raise NotImplementedError
        return Calculator.get_property(self, 'hirsh_volrat', self.atoms)
