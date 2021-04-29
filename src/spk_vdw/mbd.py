import os
import numpy as np
from ase import atoms
import ase.units as units
from ase.calculators.calculator import Calculator, FileIOCalculator
try:
    from pymbd import mbd_energy, from_volumes
    from pymbd.fortran import MBDGeom
except:
    raise ImportError("MBD calculator relies on pymbd. Please install libmbd and pymbd!")

class MBD(FileIOCalculator):

    """ MBD calculator class.

    This is an ASE interface to libmbd (https://github.com/libmbd/libmbd),
    which implements different dispersion correction methods including vdW(TS), vdW^surf, MBD, MBD-NL
    Force, energy and stress calculations are supported.
    Written by Reinhard J. Maurer, r.maurer@warwick.ac.uk
    March 2021
    """

    name = 'MBD'
    implemented_properties = ['energy', 'forces', 'stress']
  
    def __init__(self, restart=None,
                 label=os.curdir, atoms=None, **kwargs):

        # default parameters
        default_parameters = dict(k_grid=None,
                                scheme='VDW',   #'VDWw' or 'MBD', default is VDW'
                                params='TS',    #either TS or TSSURF
                                nfreq=15,
                                beta=0.83,   #PBE default value
                                ts_sr=0.94,   #PBE default value
                                do_rpa=False,
                                )

        # pound sign designates non-defaulting parameters
        valid_args = (
                'k_grid',
                'scheme',
                'params',
                'nfreq',
                'beta',
                'ts_sr',
                'do_rpa')
        
        for arg, val in default_parameters.items():
            setattr(self, arg, val)

        # set any additional keyword arguments and overwrite defaults
        # for arg, val in self.parameters.iteritems():
        for arg, val in kwargs.items():
            if arg in valid_args:
                if arg == 'params' or arg == 'scheme':
                    val = val.upper()
                    if arg == 'scheme' and val == 'TSSURF':
                        val = 'TSsurf'
                setattr(self, arg, val)
            else:
                raise RuntimeError('unknown keyword arg "%s" : not in %s'
                                   % (arg, valid_args))

        self.hirshvolrat_is_set = False
        self.hirsh_volrat = None
        
        Calculator.__init__(self, restart,
                            label, atoms, **kwargs)

    def calculate(self, atoms, properties=["energy"], 
            system_changes):
        """ actual calculation of all properties. """



        self.atoms = atoms.copy()
    
        if all(atoms.get_pbc()):
            lattice = np.array(atoms.get_cell()) / units.Bohr
        else:
            lattice = None

        if not self.hirshvolrat_is_set:
            self.hirsh_volrat = np.ones(len(atoms),dtype=np.float)

        if self.calculation_required(atoms, properties):
            Calculator.calculate(self, atoms)

            self.alpha_0, self.C6, self.R_vdw = from_volumes(
                atoms.get_chemical_symbols(), 
                self.hirsh_volrat, 
                kind=self.params
                )

            if 'forces' in properties:
                do_force = True
            else:
                do_force = False

            mbdgeom = MBDGeom(
                    coords=atoms.positions / units.Bohr, 
                    lattice=lattice, 
                    k_grid=self.k_grid,
                    n_freq=self.nfreq,
                    do_rpa=self.do_rpa,
                    )

            if self.scheme == 'MBD':
                energy = mbdgeom.mbd_energy(
                    self.alpha_0, 
                    self.C6, 
                    self.R_vdw,
                    beta=self.beta,
                    force=do_force
                    )
            elif self.scheme == 'VDW':
                energy = mbdgeom.ts_energy(
                    self.alpha_0, 
                    self.C6, 
                    self.R_vdw,
                    sR=self.ts_sr,
                    d=20.0,
                    force=do_force
                    ) 
            else: 
                raise ValueError("mbd: scheme needs to be MBD or VDW")

            self.results['energy'] = energy[0] * units.Hartree

            if do_force:
                gradients = energy[1] * units.Hartree / units.Bohr
                self.results['forces'] = -gradients
                if 'stress' in properties and all(self.atoms.get_pbc()):
                    lattgradients = energy[2] * units.Hartree / units.Bohr
                    stress = np.dot(atoms.get_cell(), lattgradients.transpose(),)+\
                        np.dot(atoms.get_positions().transpose(),gradients)
                    stress = stress / (atoms.get_volume())
                    self.results['stress'] = stress

    def set_hirshfeld(self, hirsh_volrat):
        self.hirshvolrat_is_set = True
        self.hirsh_volrat = hirsh_volrat
