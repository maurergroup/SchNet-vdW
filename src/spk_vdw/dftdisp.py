from . import sdc
#import switch
import os
import numpy as np
from ase import atoms
from ase.calculators.calculator import Calculator


class dftdisp(Calculator):

    """ DFT_DISP calculator class.

    This is an ASE interface to the semp_disc_corr.f90 dispersion
    correction implementation for DFT which includes the common
    Grimme-Correction as well as TS and TSSurf-models. Only force,
    energy and stress calculations are supported. This interface was
    written by Georg Michelitsch georg.michelitsch@theo.ch.tum.de at
    TU Munich in Dec 2013.

    """

    implemented_properties = ['energy', 'forces', 'stress']

    # default parameters
    default_parameters = dict(logfile_name='sdc_recode.log                                    ',
                              sedc_xc='PBE     ',
                              sedc_scheme='TS      ',
                              sedc_print_level=3,
                              sedc_do_pbc=-1,
                              sedc_do_num_f=0,
                              sedc_do_num_stress=0,
                              sedc_pbc_switch='ABC',
                              sedc_pbc_energy_tol=1e-6,
                              sedc_pbc_force_tol=1e-7,
                              sedc_pbc_img_fixed_nshells=0,
                              sedc_pbc_n_shells=0,
                              sedc_pbc_backfold_coord=0,
                              sedc_pbc_file_read=-1,
                              sedc_do_standalone=False)

    # pound sign designates non-defaulting parameters
    # S.R. = strictly required parameters for calculation
    valid_args = ('logfile_name',
                  'internal_cart_coord',            #
                  'sedc_cell_vectors',              #
                  'sedc_species',                   #
                  'sedc_n_ions',                    #
                  'sedc_cart_coord',                #
                  'sedc_scheme',
                  'sedc_xc',
                  'sedc_print_level',
                  'sedc_do_pbc',
                  'sedc_do_num_f',
                  'sedc_do_num_stress',
                  'sedc_pbc_switch',
                  'sedc_pbc_energy_tol',
                  'sedc_pbc_force_tol',
                  'sedc_pbc_img_fixed_nshells',
                  'sedc_pbc_n_shells',
                  'sedc_pbc_backfold_coord',
                  'sedc_pbc_g_fold',                #
                  'sedc_n_groups',                  # S.R.
                  'sedc_groups',                    # S.R.
                  'sedc_pbc_g_switches',            # S.R.
                  'sedc_pbc_g_only_intra',          #
                  'sedc_pbc_g_cells',               #
                  'sedc_ts_veff_div_vfree',         #
                  'sedc_tssurf_vfree_div_vbulk',
                  'sedc_skip_atom',                 #
                  'sedc_pbc_g_skip',                #
                  'sedc_pbc_file_read',
                  'sedc_do_standalone')

    sdch = 'sdc.sdc_recode.'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=os.curdir, atoms=None, **kwargs):

        # set any additional keyword arguments
        # for arg, val in self.parameters.iteritems():
        for arg, val in kwargs.items():
            if arg in self.valid_args:
                # satisfy our static 8-char requirement in fortran
                if (arg == ('sedc_scheme' or 'sedc_xc')):
                    setattr(self, arg, str(val).ljust(8))
                elif (arg == 'logfile_name'):
                    setattr(self, arg, str(val).ljust(50))
                else:
                    setattr(self, arg, val)
            else:
                raise RuntimeError('unknown keyword arg "%s" : not in %s'
                                   % (arg, self.valid_args))

        self.hirshvolrat_is_set = False

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)

    def get_potential_energy(self, atoms=None):
        self.update_properties(atoms)
        return eval(self.sdch + 'sedc_energy')

    def get_forces(self, atoms=None):
        self.update_properties(atoms)

        forces = eval(self.sdch + 'sedc_forces.transpose()')
        ###REMOVING TRANSLATIONS
        unit_vec = np.zeros([len(forces)*3,3])
        unit_vec[:, 0] = 1.0
        unit_vec[:,0] /= np.linalg.norm(unit_vec[:,0])
        forces -= (unit_vec[:,0] * np.dot(forces.flatten(),unit_vec[:,0])).reshape(-1,3)
        unit_vec = np.zeros([len(forces)*3,3])
        unit_vec[:, 1] = 1.0
        unit_vec[:,1] /= np.linalg.norm(unit_vec[:,1])
        forces -= (unit_vec[:,1] * np.dot(forces.flatten(),unit_vec[:,1])).reshape(-1,3)
        unit_vec = np.zeros([len(forces)*3,3])
        unit_vec[:, 2] = 1.0
        unit_vec[:,2] /= np.linalg.norm(unit_vec[:,2])
        forces -= (unit_vec[:,2] * np.dot(forces.flatten(),unit_vec[:,2])).reshape(-1,3)

        return forces

    def get_stress(self, atoms=None):
        self.update_properties(atoms)
        return eval(self.sdch + 'sedc_stress')

    def update_properties(self, atoms):
        """ check if already computed everything for this set of atoms. """

        if not hasattr(self, 'atoms') or self.atoms != atoms:
            self.calculate(atoms)

    def calculate(self, atoms):
        """ actual calculation of all properties. """

        self.atoms = atoms.copy()

        self.initialize_sdc(atoms=self.atoms)

        # Run the calculation
        sdc.sdc_recode.sedc()

    def initialize_sdc(self, atoms=None):
        """ Initialization of all parameters in the sedc-module.

        Here we initialize all parameters in order to successfully calculate
        with sdc_recode module. The priority of parameters is user > default >
        dummy. The bare minimum of required parameters is:

        * atoms             [ ASE atoms-object of the system treated        ]
        * sedc_n_groups     [ integer, number of differently treated groups ]
        * sedc_groups       [ array of integers, number of atoms per group  ]
        * sedc_pbc_switches [ np.array of 6-dim vectors specifying the VdW
                              contributions in A,B,C in positive and
                              negative direction                            ]

        This has to be done after the creation of our calculator-object,
        because otherwise we do not have the atoms-module to calculate most of
        the required properties.

        """
        for arg in self.valid_args:
            if hasattr(self, arg):
                # In order to avoid lengthy np-initialization in aims-script
                if (arg == 'sedc_pbc_g_switches'):
                    setattr(sdc.sdc_recode, arg,
                            eval('np.transpose(self.' + arg + ')'))
                elif (arg == 'sedc_tssurf_vfree_div_vbulk'):
                    tmp_array = \
                        np.array(([1] * atoms.get_number_of_atoms()),
                                np.float64)
                    species = []
                    n_species = 0
                    for atom in atoms:
                        species_i = atom.number
                        if species_i not in species:
                            species.append(species_i)
                            n_species += 1
                    species.sort()
                    for i,a in enumerate(atoms):
                        for ss, s in enumerate(species):
                            if s==a.number:
                                tmp_array[i] = \
                                    self.sedc_tssurf_vfree_div_vbulk[ss]
                    sdc.sdc_recode.sedc_tssurf_vfree_div_vbulk = \
                            tmp_array
                else:
                    setattr(sdc.sdc_recode, arg, eval('self.' + arg))
            elif arg in self.default_parameters:
                setattr(sdc.sdc_recode, arg, self.default_parameters[arg])
            else:
                if (arg == 'internal_cart_coord'):
                    sdc.sdc_recode.internal_cart_coord = \
                        atoms.get_positions().transpose().copy()
                elif (arg == 'sedc_cart_coord'):
                    sdc.sdc_recode.sedc_cart_coord = \
                        atoms.get_positions().transpose().copy()
                elif (arg == 'sedc_cell_vectors'):
                    sdc.sdc_recode.sedc_cell_vectors = \
                        atoms.get_cell().transpose().copy()
                elif (arg == 'sedc_species'):
                    sdc.sdc_recode.sedc_species = \
                        atoms.get_atomic_numbers().copy()
                elif (arg == 'sedc_n_ions'):
                    sdc.sdc_recode.sedc_n_ions = \
                        atoms.get_number_of_atoms()
                elif (arg == 'sedc_pbc_g_fold'):
                    sdc.sdc_recode.sedc_pbc_g_fold = [0] * self.sedc_n_groups
                elif (arg == 'sedc_ts_veff_div_vfree'):
                    if not self.hirshvolrat_is_set:
                        sdc.sdc_recode.sedc_ts_veff_div_vfree = \
                            np.array(([1] * atoms.get_number_of_atoms()),
                                     np.float64)
                elif (arg == 'sedc_tssurf_vfree_div_vbulk'):
                    sdc.sdc_recode.sedc_tssurf_vfree_div_vbulk = \
                        np.array(([1] * atoms.get_number_of_atoms()),
                                np.float64)
                elif (arg == 'sedc_skip_atom'):
                    sdc.sdc_recode.sedc_skip_atom = \
                        np.array(([-1] * atoms.get_number_of_atoms()),
                                 np.float64)
                elif (arg == 'sedc_pbc_g_only_intra'):
                    sdc.sdc_recode.sedc_pbc_g_only_intra = \
                        [0] * self.sedc_n_groups
                elif (arg == 'sedc_pbc_g_cells'):
                    sdc.sdc_recode.sedc_pbc_g_cells = \
                        np.tile(atoms.get_cell().transpose().copy(),
                                (1, self.sedc_n_groups))
                elif (arg == 'sedc_pbc_g_skip'):
                    sdc.sdc_recode.sedc_pbc_g_skip = [0] * self.sedc_n_groups
                else:
                    print("You've been sloppy my friend. :) Variable:", arg, \
                          " does not exist!")


    # TODO: [ ] Is this the correct value to write to?
    def set_hirshfeld(self, hirsh_volrat):
        self.hirshvolrat_is_set = True
        setattr(sdc.sdc_recode, 'sedc_ts_veff_div_vfree', hirsh_volrat)
