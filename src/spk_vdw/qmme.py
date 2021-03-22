from ase.calculators.calculator import Calculator
from ase import Atoms
from ase.constraints import FixAtoms
import numpy as np
import os
from time import strftime, gmtime


class qmme(Calculator):
    """ QM-Me calculator.

    This is a very generic QM/MM code for ASE, based
    upon the very specialized qmmm_manyqm code present in ASE.

    """

    implemented_properties = ['energy', 'forces']

    valid_args = ('nqm_regions',        # int
                  'nmm_regions',        # int
                  'qm_calculators',     # array of calculators
                  'qm_atoms',           # array of tuples
                  'mm_calculators',     # array of calculators
                  'mm_atoms',           # array of tuples
                  'mm_mode',            # string
                  'qm_pbcs',            # array of 3x1 boolean arrays
                  'mm_pbcs',            # array or 3x1 boolean arrays
                  'qm_cell',            # array or 3x3 arrays
                  'mm_cell',            # array of 3x3 arrays
                  'hirbulk',            # array of floats
                  'hirlast',            # array of tuples
                  'freeze')             #freeze atoms #Changed JW: added

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=os.curdir, atoms=None, logprfx=None,
                 reuse=False, hirlog=False, reset=False, **kwargs):
        #Changed JW. set reset to False as spk-calculators do not have reset
        """ Assign QM and MM calculators and regions.

        Parameters
        ==========
        nqm_regions:        int
            how many qm regions

        nmm_regions:        int
            how many MM regions, defaults to 1.

        qm_calculators:     list of members of a Class defining a Calculator
            ase qm-calculator for each qm region

        mm_calculators:      member of a Class defining a Calculator
            ase mm-calculator for the mm region (the whole system)

        qm_atoms:           array of tuples
            areas for the QM-calculation

        mm_atoms:           array of tuples
            areas for the MM-calculation

        mm_mode:            string ('allatoms', 'complementary, 'explicit')
            specifies how the mm_region will be initialized

        qm_pbcs:             array or 3-dim boolean vectors
            periodic boundary conditions for the QM regions

        mm_pbcs:             array or 3-dim boolean vectors
            periodic boundary conditions for the MM regions

        qm_cell:            array of 3x3 arrays
            unit cell parameters for each QM region

        mm_cell:            array of 3x3 arrays
            unit cell parameters for each MM region

        hirbulk:            list of floats
            hirshfeld ratio value for each of nhf_regions

        hirlast:            list of 2-tuples
            first and last atom belonging to each in nhf_regions

        """

        # Set any keyword arguments
        for arg, val in kwargs.items():
            if arg in self.valid_args:
                setattr(self, arg, val)
            else:
                raise RuntimeError('unknown keyword arg "%s" : not in %s'
                                   % (arg, self.valid_args))

        # Check the user input for missing information
        error_head = ' +----------------** QMME WARNING **----------------+'
        error_tail = ' +--------------------------------------------------+'
        for arg in self.valid_args:
            if arg not in kwargs.keys():
                if arg == 'nqm_regions':
                    if hasattr(self, 'qm_atoms'):
                        print(error_head)
                        print(' |  Keyword  nqm_regions  not specified,')
                        print(' |  defaulting to len(qm_atoms).')
                        print(error_tail + '\n')
                        self.nqm_regions = len(self.qm_atoms)
                    else:
                        print(error_head)
                        print(' |  Keyword  nqm_regions  not specified,')
                        print(' |  defaulting to 0. No QM-regions detected!')
                        print(error_tail + '\n')
                        self.nqm_regions = 0
                if arg == 'nmm_regions':
                    if hasattr(self, 'mm_atoms'):
                        print(error_head)
                        print(' |  Keyword  nmm_regions  not specified,')
                        print(' |  defaulting to len(mm_atoms).')
                        print(error_tail + '\n')
                        self.nmm_regions = len(self.mm_atoms)
                    elif hasattr(self, 'mm_mode') and \
                            hasattr(self, 'mm_calculators'):
                            print(error_head)
                            print(' |  Keyword  nmm_regions  not specified,')
                            print(' |  defaulting to len(mm_calculators).')
                            print(error_tail + '\n')
                            self.nmm_regions = len(self.mm_calculators)
                    else:
                        print(error_head)
                        print(' |  Keyword  nmm_regions  not specified,')
                        print(' |  defaulting to 0. No MM-regions detected!')
                        print(error_tail + '\n')
                        self.nmm_regions = 0
                if arg == 'qm_calculators':
                    if hasattr(self, 'qm_atoms'):
                        raise RuntimeError(
                            'WARNING: Keyword  qm_calculators  not ' +
                            'specified! Please specify calculators ' +
                            'for QM regions.\n')
                    else:
                        break
                if arg == 'mm_calculators':
                    if hasattr(self, 'mm_atoms') or hasattr(self, 'mm_mode'):
                        raise RuntimeError(
                            'WARNING: Keyword  mm_calculators  not ' +
                            'specified! Please specify calculators ' +
                            'for MM regions.\n')
                    else:
                        break
                if arg == 'qm_atoms':
                    if self.nqm_regions != 0 or \
                       hasattr(self, 'qm_calculators'):
                        raise RuntimeError(
                            'WARNING: Keyword  qm_atoms  not ' +
                            'specified! Please specify atoms ' +
                            'for QM regions.\n')
                if arg == 'mm_atoms':
                    if not hasattr(self, 'mm_mode') and self.nmm_regions != 0:
                        raise RuntimeError(
                            'WARNING: Keyword  mm_atoms  not ' +
                            'specified! Please specify atoms ' +
                            'for MM regions.\n')
                if arg == 'mm_mode':
                    if hasattr(self, 'mm_calculators') or \
                            self.nmm_regions != 0:
                        raise RuntimeError(
                            'WARNING: Keyword  mm_mode  not' +
                            'specified! Please specify mm_mode as either' +
                            'allatoms, complementary or explicit.\n')
                if arg == 'qm_pbcs':
                    print(error_head)
                    print(' |  Keyword  qm_pbcs  not specified, QM regions will')
                    print(' |  have the periodic boundary conditions of the')
                    print(' |  whole system assigned to them.')
                    print(error_tail + '\n')
                if arg == 'mm_pbcs':
                    print(error_head)
                    print(' |  Keyword  mm_pbcs  not specified, MM regions will')
                    print(' |  have the periodic boundary conditions of the')
                    print(' |  whole system assigned to them.')
                    print(error_tail + '\n')
                if arg == 'qm_cell':
                    print(error_head)
                    print(' |  Keyword  qm_cell  not specified, QM regions will')
                    print(' |  have the unit cell parameters of the')
                    print(' |  whole system assigned to them.')
                    print(error_tail + '\n')
                if arg == 'mm_cell':
                    print(error_head)
                    print(' |  Keyword  mm_cell  not specified, MM regions will')
                    print(' |  have the unit cell parameters of the')
                    print(' |  whole system assigned to them.')
                    print(error_tail + '\n')

        # Force user to specify atoms-object at QMMM-initialization
        # At the same time deliberately disable constraints within
        # our atoms.object as Reinhard Maurer suggested.
        try:
            atoms.positions
        except:
            raise RuntimeError('Keyword  atoms  not specified! QMMM-module ' +
                               'strictly requires the atoms ' +
                               'object at initialization.\n')

        if logprfx is None and hirlog:
            print(error_head)
            print(' |  Keyword  logprfx  not specified! We will create')
            print(' |  uniquely named logdirs for every calculator')
            print(' |  and each iteration step.')
            print(' |  Warning: This might make it difficult')
            print(' |           to match logfiles to calculations.')
            print(error_tail + '\n')

            self.logprfx = strftime("%b-%d_%H%M%S_", gmtime()) + 'QMMM'
        else:
            self.logprfx = logprfx

        # Routines for logging and reuse of prior data
        self.reuse = reuse

        # reset to v_hirsh/v_free = 1 if not give
        if not hasattr(self, 'hirbulk'):
            self.hirbulk = [1.0]

        if not self.reuse and hirlog:
            try:
                os.makedirs(self.logprfx)

                for q in range(self.nqm_regions):
                    os.makedirs(self.logprfx + '/QM' + str(q) + '_' +
                                self.qm_calculators[q].__class__.__name__)

                for m in range(self.nmm_regions):
                    os.makedirs(self.logprfx + '/MM' + str(m) + '_' +
                                self.mm_calculators[m].__class__.__name__)
            except OSError:
                raise RuntimeError(
                    '\nWARNING: Desired log-directories are already present:\n' +
                    'Please remove or backup old logs before restarting a \n' +
                    'calculation. Alternatively choose a different logprfx. \n' +
                    'If you wish to reuse data from previous calculations, \n'
                    'please set  reuse=True  as an additional keyword. \n')

        # Check if the QM-calculator support Hirshfeld-paritioning and inform
        # user about missing Hirshfeld-Treatment.
        if 'dftdisp' in str(self.mm_calculators):
            for qmcalc in self.qm_calculators:
                 #if "read_hirsh_volrat" not in dir(qmcalc):
                if not hasattr(qmcalc, 'get_hirsh_volrat'):
                    print(error_head)
                    print(' |  QM calculator ' + str(qmcalc.__class__.__name__) + \
                          ' does not support Hirshfeld-')
                    print(' |  partitioning. Defaulting to v_hirsh/v_free = bulk_value')
                    print(error_tail + '\n')

        # Print a Summary of the Configuration of QMMM-Calculator
        print (' +----------------** QMME SUMMARY **----------------+\n |')
        if hirlog:
            print( ' |  QMMM calculator  ' + self.logprfx + '  initialized:')
        for n in range(self.nqm_regions):
            print( ' |    QM region %i: %s (atoms: %s)' \
                % (n, self.qm_calculators[n].__class__.__name__,
                    str(self.qm_atoms[n]).strip('[]')))
            try:
                print( ' |       PBCs     : %s' \
                    % (str(self.qm_pbcs[n]).strip('[]')))
            except:
                print( ' |       PBCs     : automatic ')
            try:
                tmp_cell = ''
                for i in range(3):
                    tmp_cell += '\n |           ' + \
                        str(["%.3f" % elem for elem in self.qm_cell[n][i]])
                print( ' |       Unit Cell: ' + tmp_cell + '\n |')
            except:
                print( ' |       Unit cell: automatic \n |')
        for n in range(self.nmm_regions):
            if self.mm_mode == 'explicit':
                print( ' |    MM region %i: %s (atoms: %s)' \
                    % (n, self.mm_calculators[n].__class__.__name__,
                        str(self.mm_atoms[n]).strip('[]')))
                try:
                    print( ' |       PBCs     : %s' \
                        % (str(self.mm_pbcs[n]).strip('[]')))
                except:
                    print( ' |       PBCs     : automatic ')
                try:
                    tmp_cell = ''
                    for i in range(3):
                        tmp_cell += '\n|           ' + \
                            str(["%.3f" % elem for elem in self.mm_cell[n][i]])
                    print( ' |       Unit Cell: ' + tmp_cell + '\n |')
                except:
                    print( ' |       Unit cell: automatic \n |')
            else:
                print( ' |    MM region %i: %s (atoms: %s)' \
                    % (n, self.mm_calculators[n].__class__.__name__,
                        self.mm_mode))
                try:
                    print( ' |       PBCs     : %s' \
                        % (str(self.mm_pbcs[n]).strip('[]')))
                except:
                    print( ' |       PBCs     : automatic ')
                try:
                    tmp_cell = ''
                    for i in range(3):
                        tmp_cell += '\n|           ' + \
                            str(["%.3f" % elem for elem in self.mm_cell[n][i]])
                    print( ' |       Unit Cell: ' + tmp_cell + '\n |')
                except:
                    print( ' |       Unit cell: automatic \n |')
        print( ' +--------------------------------------------------+\n')

        # Initialization of empty data structures of calculation results.
        # This is required for some internal reconstruction purposes and
        # error processing.
        self.qm_energies = [[0.0] for i in range(self.nqm_regions)]
        self.qm_forces = [[[0.0] * atoms.get_number_of_atoms()] * 3
                          for i in range(self.nqm_regions)]
        self.mm_energies = [[0.0] for i in range(self.nmm_regions)]
        self.mm_forces = [[[0.0] * atoms.get_number_of_atoms()] * 3
                          for i in range(self.nmm_regions)]
        self.energy = 0.0
        self.forces = np.array(([[0.0] * 3] *
                                atoms.get_number_of_atoms()))
        self.qm_regions = [None for i in range(self.nqm_regions)]
        self.mm_regions = [None for i in range(self.nmm_regions)]
        self.qm_hirshpart = [0.0 for i in range(atoms.get_number_of_atoms())]

        self.nopt = 0
        self.hirlog = hirlog
        # Set the reset-value
        self.reset = reset

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms=atoms, **kwargs)

    def get_potential_energy(self, atoms=None,force_consistent=True):
        self.update_properties(atoms)
        return self.energy

    def get_forces(self, atoms=None):
        self.update_properties(atoms)
        return self.forces

    def update_properties(self, atoms):
        """ check if already computed everything for this set of atoms. """

        # if not hasattr(self, 'atoms') or self.atoms != atoms:
        # just calculate it - have the calculators below check if they already
        # know these coordinates, etc
        self.generate_calculators(atoms)
        self.calculate()

    def generate_calculators(self, atoms):

        self.atoms = atoms.copy()

        # Generate the QM/MM internal maps
        self.generate_qmmm_maps(atoms=self.atoms)

        # Initialize the atoms-objects for each QM and MM calculator
        if self.check_if_set_and_value('nqm_regions') and \
           self.check_if_set_and_value('qm_atoms'):
            for iregion in range(len(self.qm_atoms)):
                #self.generate_qm_region(region=iregion)
                # Append this new created region to the qm_regions - array
                self.qm_regions[iregion] = \
                    self.generate_qm_region(region=iregion).copy()

                # If specified, set the cells and PBCs for each region,
                # otherwise default to the initial slabs properties
                try:
                    self.qm_regions[iregion].set_pbc(self.qm_pbcs[iregion])
                except:
                    self.qm_regions[iregion].set_pbc(self.atoms.get_pbc())
                try:
                    self.qm_regions[iregion].set_cell(self.qm_cell[iregion])
                except:
                    self.qm_regions[iregion].set_cell(self.atoms.get_cell())

                # Set the calculator for this region
                self.qm_regions[iregion].set_calculator(
                    self.qm_calculators[iregion])

        # Initialize the atoms-objects of the MM regions
        if self.check_if_set_and_value('nmm_regions'):
            for iregion in range(self.nmm_regions):
                self.mm_regions[iregion] = \
                    self.generate_mm_region(region=iregion).copy()

                # If specified, set the cells and PBCs for each region,
                # otherwise default to the initial slabs properties
                try:
                    self.mm_regions[iregion].set_pbc(self.mm_pbcs[iregion])
                except:
                    self.mm_regions[iregion].set_pbc(self.atoms.get_pbc())
                try:
                    self.mm_regions[iregion].set_cell(self.mm_cell[iregion])
                except:
                    self.mm_regions[iregion].set_cell(self.atoms.get_cell())

                # Set the calculator for this region
                self.mm_regions[iregion].set_calculator(
                    self.mm_calculators[iregion])

    def calculate(self):
        """ actual calculation of all properties. """

        owd = os.getcwd()

        # Run the calculation
        if self.check_if_set_and_value('nqm_regions') and \
           self.check_if_set_and_value('qm_atoms'):
            for iqmreg, qmreg in enumerate(self.qm_regions):

                # reset the calculator
                if self.reset:
                    self.qm_calculators[iqmreg].reset()

                # Jump in the appropriate directory and start calculation
                if self.hirlog:
                    try:
                        curdir = self.logprfx + '/QM' + str(iqmreg) + '_' + \
                            self.qm_calculators[iqmreg].__class__.__name__
                        os.chdir(curdir)
                    except OSError:
                        os.makedirs(curdir)
                        os.chdir(curdir)

                    if not self.reuse or \
                       not os.path.exists('run' + str(self.nopt)):
                        os.makedirs('run' + str(self.nopt))
                    os.chdir('run' + str(self.nopt))

                qmreg.set_calculator(self.qm_calculators[iqmreg])
                self.qm_energies[iqmreg] = qmreg.get_potential_energy()
                self.qm_forces[iqmreg] = qmreg.get_forces()

                # Jump back to the root directory
                os.chdir(owd)

        if self.check_if_set_and_value('nmm_regions'):
            for immreg, mmreg in enumerate(self.mm_regions):

                # reset the calculator
                if self.reset:
                    self.mm_calculators[immreg].reset()

                # Jump in the appropriate directory and start calculation
                if self.hirlog:
                    try:
                        curdir = self.logprfx + '/MM' + str(immreg) + '_' + \
                            self.mm_calculators[immreg].__class__.__name__
                        os.chdir(curdir)
                    except OSError:
                        os.makedirs(curdir)
                        os.chdir(curdir)

                    if not self.reuse or \
                       not os.path.exists('run' + str(self.nopt)):
                        os.makedirs('run' + str(self.nopt))
                    os.chdir('run' + str(self.nopt))

                try:
                    self.mm_calculators[immreg].set_hirshfeld(
                        self.generate_hirshpart(iregion=immreg) )
                except AttributeError:
                    pass
                    ## Be silent.
                    #print self.mm_calculators[immreg].__class__.__name__ + \
                        #"calculator does not yet support Hirshfeld -\n" + \
                        #"partitioning assisted calculation. This is only\n" + \
                        #"a warning about this fact. The calculation\n" + \
                        #"will proceed uninterruptedly.\n"
                mmreg.set_calculator(self.mm_calculators[immreg])
                self.mm_energies[immreg] = mmreg.get_potential_energy()
                self.mm_forces[immreg] = mmreg.get_forces()

                # Jump back to the root directory
                os.chdir(owd)

        self.nopt += 1

        self.energy = sum(self.qm_energies) + sum(self.mm_energies)

        # Since our force vector was split into parts, reconstruct it here
        self.reconstruct_force_vector()

    def generate_qm_region(self, region):
        """ Generates QM-regions and initializes calculators.

        Here we generate the qm_regions objects and pass them to the selected
        calculator. The region input is expected to be passed as an array of
        tuples assigning subregions of the whole qm_region per calculator,
        even allowing non-continguous assignment.

        Example
        =======
        qm_atoms = ([(2,7)],[(0,2),(4,6)])
                This defines two qm_regions with atoms 2-7 in the first region
                and 0-2 and 4-6 in the second region. Each region will be
                assigned to the corresponding calculator in qm_calculators.

        """

        # Initialize empty atoms object and array
        qm_region = Atoms() 

        # Assign selected slices to new atoms object
        for subregion in range(len(self.qm_atoms[region])):
            #changed JW does not work otherwise
            if subregion == 0:
                qm_region = self.atoms[slice(self.qm_atoms[region][subregion][0],self.qm_atoms[region][subregion][1])]
            else:
                qm_region += self.atoms[slice(self.qm_atoms[region][subregion][0],self.qm_atoms[region][subregion][1])]
        
        # Return the newly created qm_region
        qm_region.set_constraint(FixAtoms(self.freeze)) #changed JW
        return qm_region
        # Append this new created region to the qm_regions - array
        #self.qm_regions.append(qm_region.copy())
        # Dec 28 2013: because of faulty variable scope this was moved to
        # calculate()
        # mit slice() slices aus der qm_atoms variable erzeugen und mit
        # atoms.extend() methode dynamisch zusammenfuddeln -> ermoeglicht
        # non-continguous areas fuer qm_regions
        # qm_region = self.atoms[self.qm_atoms[region]]

    # TODO: [ ] Error handling, if certain variables are not set
    def generate_mm_region(self, region):
        """ Generates MM-regions and initializes calculators.

        This basically does the same thing like the corresponding qm-generator
        but it allows different ways of initialization depending on the status
        of self.mm_mode .

        Parameters
        ==========
        mm_mode:    string, either 'allatoms', 'complementary' or 'explicit'

            allatoms:       MM calculator is initialized with all atoms
                            of the initial atoms-object, default
                            Recommended for systems with only one MM calculator

            complementary:  MM calculator is initialized with all atoms but
                            those declared as 'qm_atoms'

            explicit:       MM calculator is initialized with only those atoms
                            specified in 'mm_atoms'
                            Required if more than one MM calculator are used
        """

        if (self.mm_mode == 'allatoms'):

        # Attach the whole atoms-object
    #            for region in range(self.nmm_regions):
            # TODO: [ ] These are very dirty defaults!! HARDCODED TODO TODO
            # TODO
            # self.mm_calculators[region].set_atoms(self.atoms)
            #mm_region.set_pbc(self.atoms.get_pbc())
            #mm_region.set_cell(self.atoms.get_cell())
            # TODO: [ ] Jan 14': This is not needed? Not sure why it was left
            # here.
            #self.atoms.copy().set_calculator(self.mm_calculators[region])

            # Append this new created region to the mm_regions - array
            #self.mm_regions.append(self.atoms.copy())
            mm_region = self.atoms.copy()
            mm_region.set_constraint(FixAtoms(self.freeze)) #changed JW
            return mm_region

        elif (self.mm_mode == 'complementary') or (self.mm_mode == 'explicit'):

            # Initialize empty atoms object
            mm_region = Atoms()

            # TODO: [x] not implemented yet. This is not a priority.
            # Basic idea: Overlay range(self.atoms) with a map of the
            #             qm_slices and create slices from the the
            #             remaining positions.

            # Construct new atoms object
            #for region in range(self.nmm_regions):
            for atom in (range(len(self.atoms))):
                #Changed JW
                if (self.mm_map[region][atom]):
                    mm_region.append(self.atoms[atom])

        # TODO: [x] Removeable code, replaced by more generic mapping-method
        #elif (self.mm_mode == 'explicit'):

            ## Initialize empty atoms object
            #mm_region = Atoms()

            ## Assign selected slices to new atoms object
            #for region in range(self.nmm_regions):
                #for subregion in range(len(self.mm_atoms[region])):
                    #mm_region +=\
                        #self.atoms[slice(
                            #self.mm_atoms[region][subregion][0],
                            #self.mm_atoms[region][subregion][1])]

            # TODO: [ ] die ham alle kane set_atoms()
            # Attach the new atoms-object (complementary or explicit)
            # TODO: [ ] Jan 14': This is not needed? Not sure why it was left
            # here.
            # self.mm_calculators[region].set_atoms(mm_region)

            # Append this newly created region to the mm_regions - array
            # self.mm_regions.append(mm_region.copy())
            mm_region.set_constraint(FixAtoms(self.freeze)) #changed JW
            return mm_region

    def generate_qmmm_maps(self, atoms):
        """ This function generates boolean maps of QM and MM regions. """

        if self.check_if_set_and_value('nqm_regions') and \
           self.check_if_set_and_value('qm_atoms'):
            self.qm_map = np.array([[False] *
                                    self.atoms.get_number_of_atoms()] *
                                   self.nqm_regions)
            for region in range(len(self.qm_atoms)):
                for subregion in range(len(self.qm_atoms[region])):
                    self.qm_map[region][slice(
                        self.qm_atoms[region][subregion][0],
                        self.qm_atoms[region][subregion][1])] = True

        if self.check_if_set_and_value('nmm_regions'):
            if (self.mm_mode == 'allatoms'):
                self.mm_map = np.array([[True] *
                                       self.atoms.get_number_of_atoms()] *
                                       self.nmm_regions)

            elif (self.mm_mode == 'complementary'):
                # This will fail if we didnt generate a qm_map. However,
                # running a QMMM calculation without QM is not logical.
                # [as in Vulcan logic]
                self.mm_map = np.array([np.invert(np.any(
                    self.qm_map.copy(), axis=0))])

            elif (self.mm_mode == 'explicit'):
                self.mm_map = np.array([[False] *
                                       self.atoms.get_number_of_atoms()] *
                                       self.nmm_regions)
                for region in range(len(self.mm_atoms)):
                    for subregion in range(len(self.mm_atoms[region])):
                        self.mm_map[region][slice(
                            self.mm_atoms[region][subregion][0],
                            self.mm_atoms[region][subregion][1])] = True

    def reconstruct_force_vector(self):
        """ This function reconstructs the force vector.

        This function generates boolean arrays defining the QM and
        MM-regions in order to be able to reconstruct the initial
        order of atoms in the atoms-object and be able to sum up
        those necessary.
        """

        self.forces = np.array(([[0.0] * 3] *
                                self.atoms.get_number_of_atoms()))
        # Reconstruct Forces with QM-Results
        if self.check_if_set_and_value('nqm_regions') and \
           self.check_if_set_and_value('qm_atoms'):
            for region in range(self.nqm_regions):
                nelem = np.sum(self.qm_map[region])
                for atom in range(len(self.atoms)):
                    if (self.qm_map[region][atom]) and (nelem > 0):
                        self.forces[atom] += self.qm_forces[region][
                            np.sum(self.qm_map[region]) - nelem]
                        nelem -= 1

        # Reconstruct Forces with MM-Results
        if self.check_if_set_and_value('nmm_regions'):
            for region in range(self.nmm_regions):
                nelem = np.sum(self.mm_map[region])
                for atom in range(len(self.atoms)):
                    if (self.mm_map[region][atom]) and (nelem > 0):
                        self.forces[atom] += self.mm_forces[region][
                            np.sum(self.mm_map[region]) - nelem]
                        nelem -= 1

    def get_hirsh_volrat(self, i_current_qm):
        """ This function returns the hirshfeld partitioning ratio.

        This function returns the partitioning ratio and relies on the
        calculator-object to be able to compute it. If this is not implemented,
        an exception will be thrown. It will return an array with floats of
        the Hirshfeld-ration v_hirsh/v_free .
        """
        if (not self.solved_hvr_qm[i_current_qm]):
            try:
                self.hvrs_qm_calc[i_current_qm] = self.qm_calculators[i_current_qm].get_hirsh_volrat()
            except AttributeError:
                # !DO NOT! be silent and take default data.
                print( self.qm_calculators[i_current_qm].__class__.__name__ + \
                    " calculator does not yet support Hirshfeld-partitioning. " + \
                    " Defaulting to v_hirsh/v_free = 1.")
                self.hvrs_qm_calc[i_current_qm] = [1.0 for i in \
                                  range(self.qm_regions[i_current_qm].get_number_of_atoms())]

            self.solved_hvr_qm[i_current_qm] = True

        return self.hvrs_qm_calc[i_current_qm]


    def generate_hirshpart(self, iregion):
        """ This function generates an array with all the hirshfeld parameters.

        Since only after the QM calculation, the Hirshfeld charges required in
        MM treatment are available for QM atoms we need to reconstruct an array
        which contains those as well as a generic value of self.hirbulk for all
        not available atoms.
        """

        if self.nqm_regions == 0:
            self.qm_hirshpart = [1.0 for i
                                 in range(self.atoms.get_number_of_atoms())]

        # Construct Hirshfeld Partitioning with QM-Results
        self.solved_hvr_qm = [False,]*self.nqm_regions
        self.hvrs_qm_calc = []
        for region in range(self.nqm_regions):
            nelem = np.sum(self.qm_map[region])
            self.hvrs_qm_calc.append(np.zeros(nelem))
            for atom in range(len(self.atoms)):
                if (self.qm_map[region][atom]) and (nelem > 0):
                    self.qm_hirshpart[atom] = self.get_hirsh_volrat(region)[
                        np.sum(self.qm_map[region]) - nelem]
                    nelem -= 1
                # set to the bulk value
                elif (self.qm_hirshpart[atom] == 0):
                    self.qm_hirshpart[atom] = self.hirbulk[0]

        # If we requested to overwrite particular atoms the
        # hirshfeld-bulk-value then is done here
        if hasattr(self, 'hirlast'):
            for nhf, hirlast in enumerate(self.hirlast):
                s = slice(hirlast[0], hirlast[1])
                self.qm_hirshpart[s] = [self.hirbulk[nhf]] * (hirlast[1] -
                                                              hirlast[0])

        ## reset calculator flag for new evaluation
        self.solved_hvr_qm = [False,]*self.nqm_regions
        # Correct Hirshfeld-Partitioning for smaller portions of the whole
        # atoms-object if explicit or complementary mode is chosen
        if self.mm_mode in ['explicit', 'complementary']:
            nelem = np.sum(self.mm_map[iregion])
            # Everthing not touched by this should actually be 0!
            qm_hirshtmp = [0.0 for i in range(np.sum(self.mm_map[iregion]))]
            for atom in (range(len(self.atoms))):
                if (self.mm_map[iregion][atom]) and (nelem > 0):
                    qm_hirshtmp[np.sum(self.mm_map[iregion]) - nelem] = \
                        self.qm_hirshpart[atom]
                    nelem -= 1
            # also add the bulk values here if requested
            if hasattr(self, 'hirlast'):
                for nhf, hirlast in enumerate(self.hirlast):
                    s = slice(hirlast[0], hirlast[1])
                    qm_hirshtmp[s] = [self.hirbulk[nhf]] * (hirlast[1] -
                                                            hirlast[0])

            return qm_hirshtmp

        # Return the final array
        return self.qm_hirshpart

    # TODO: [ ] This function can be removed if all for-loops with a
    #           range-argument over nqm_regions or nmm_regions are
    #           correctly evaluated and therefore not looped if 0.
    def check_if_set_and_value(self, string):
        if hasattr(self, string):
            if string == 'qm_atoms':
                if self.qm_atoms[0]:
                    return True
            elif (eval('self.' + string) > 0):
                return True
        return False
