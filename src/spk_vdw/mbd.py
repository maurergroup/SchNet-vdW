import os
import numpy as np
from ase import atoms
import ase.units as units
from ase.calculators.calculator import Calculator
try:
    from pymbd import mbd_energy_species, from_volumes
    from pymbd.fortran import MBDGeom
except:
    raise ImportError("MBD calculator relies on pymbd. Please install libmbd and pymbd!")

class MBD(Calculator):

    """ MBD calculator class.

    This is an ASE interface to libmbd (https://github.com/libmbd/libmbd),
    which implements different dispersion correction methods including vdW(TS), vdW^surf, MBD, MBD-NL
    Force, energy and stress calculations are supported.
    Written by Reinhard J. Maurer, r.maurer@warwick.ac.uk
    March 2021
    """

    implemented_properties = ['energy', 'forces', 'stress']

    # default parameters
    default_parameters = dict(k_grid=None,
                              scheme='VDW',   #'VDWw' or 'MBD', default is VDW'
                              params='TS',    #either TS or TSSURF
                              n_freq=15,
                              beta='0.83',   #PBE default value
                              ts_sr='0.94',   #PBE default value
                              do_rpa=False,
                              )

    # pound sign designates non-defaulting parameters
    valid_args = (
            'k_grid',
            'scheme',
            'params',
            'n_freq',
            'beta',
            'ts_sr',
            'do_rpa')

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=os.curdir, atoms=None, **kwargs):

        for arg, val in default_parameters:
            setattr(self, arg, val)

        # set any additional keyword arguments and overwrite defaults
        # for arg, val in self.parameters.iteritems():
        for arg, val in kwargs.items():
            if arg in self.valid_args:
                if arg == 'params' or arg == 'scheme':
                    arg = arg.upper()
                setattr(self, arg, val)
            else:
                raise RuntimeError('unknown keyword arg "%s" : not in %s'
                                   % (arg, self.valid_args))

        self.hirshvolrat_is_set = False
        self.hirsh_volrat = np.ones(len(atoms))
        self.update_mbd(atoms)
        
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)

    def update_properties(self, atoms):
        """ check if already computed everything for this set of atoms. """

        if not hasattr(self, 'atoms') or self.atoms != atoms:
            self.calculate(atoms)

    def calculate(self, atoms):
        """ actual calculation of all properties. """

        self.atoms = atoms.copy()
        self.update_mbd(atoms=self.atoms)
        if self.scheme == 'MBD':
            energy, force, stress = self.mbdgeom.mbd_energy(
                self.alpha_0, 
                self.C6, 
                self.R_vdw,
                beta=self.beta,
                force=True
                )
        elif self.scheme == 'VDW':
            energy, force, stress = self.mbdgeom.ts_energy(
                self.alpha_0, 
                self.C6, 
                self.R_vdw,
                sR=self.ts_sr,
                d=20.0,
                force=True
                ) 
        
       self.results['energy'] = energy * units.Hartree
       self.results['forces'] = -force * units.Hartree / units.Bohr
       self.results['stress'] = -np.dot(atoms.get_cell(),stress) * units.Hartree / units.Bohr


    def update_mbd(self, atoms=None):
        """ Initialization of libmbd geom_t object via exposed MBDGeom.

        """
        
        if all(atoms.get_pbc()):
            lattice = atoms.get_cell()
        elif any(atoms.get_pbc()):
            lattice = None # mixed pbc currently not implemented TODO
        else:
            lattice = None

        self.alpha_0, self.C6, self.R_vdw = from_volumes(atoms.get_chemical_symbols(), self.hirsh_volrat, kind=self.params)

        self.mbdgeom = MBDGeom(
                coords=atoms.positions, 
                lattice=lattice, 
                k_grid=self.k_grid,
                n_freq=self.nfreq,
                do_rpa=self.do_rpa,
                )


    def set_hirshfeld(self, hirsh_volrat):
        self.hirshvolrat_is_set = True
        self.hirsh_volrat = hirsh_volrat
