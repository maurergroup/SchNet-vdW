import os
import subprocess
import numpy as np
from ase.calculators.calculator import FileIOCalculator, Calculator


class mbdio(FileIOCalculator):

    """ DFT_MBD IO-calculator class.

    This is an simple I/O ASE interface to the MBD_AT_rsSCS.F90
    many body dispersion implementation (A. Tkatchenko PRL 108,
    236402 (2012)). The given code was obtained from R. Maurer and
    is used with permission. Forces will be evaluated numerically.
    Energy is parsed from stdout. This interface was
    written by Georg S. Michelitsch georg.michelitsch@theo.ch.tum.de at
    TU Munich in Jun 2015.

    """

    implemented_properties = ['energy', 'forces']
    changes = ['positions', 'numbers', 'cell', 'pbc']

    # default parameters
    #  xc 1   #(for PBE type)
    #  xc 2   #(for PBE0 type)
    #  xc 3   #(for HSE type)
    default_parameters = dict(INPUT_FILE='benzene.xyz',
                              SETTING_FILE='setting1.in',
                              OUTPUT_FILE='MBD.log',
                              xc='1',
                              mbd_cfdm_dip_cutoff='200.d0',
                              mbd_scs_dip_cutoff='200.0',
                              mbd_supercell_cutoff='15.d0',
                              mbd_scs_vacuum_axis='.true. .true. .false.',
                              mode='TS') #default otherwise rsSCS

    valid_args = ('INPUT_FILE',
                  'SETTING_FILE',
                  'OUTPUT_FILE',
                  'xc',
                  'mbd_cfdm_dip_cutoff',
                  'mbd_scs_dip_cutoff',
                  'mbd_supercell_cutoff',
                  'mbd_scs_vacuum_axis')

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, command=None, calculate_forces=True,
                 **kwargs):
        self.command = command
        self.calculate_forces = calculate_forces
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise RuntimeError('Please set $%s environment variable ' %
                               ('ASE_' + self.name.upper() + '_COMMAND') +
                               'or supply the command keyword')
        olddir = os.getcwd()
        try:
            syscall = self.command + ' ' + self.parameters['INPUT_FILE'] + \
                ' ' + self.parameters['SETTING_FILE'] + ' > ' + \
                self.parameters['OUTPUT_FILE']
            errorcode = subprocess.call(syscall, shell=True)
        finally:
            os.chdir(olddir)

        if errorcode:
            raise RuntimeError('%s returned an error: %d' %
                               (self.name, errorcode))
        self.read_results()

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input file(s).

        Call this method first in subclasses so that directories are
        created automatically."""

        # Write geometry-file
        b = open(self.parameters['INPUT_FILE'], 'w+')
        b.write(str(atoms.get_number_of_atoms()) + '\n\n')
        for a in range(len(atoms)):
            pos = ''
            for i in range(3):
                pos += str(atoms.get_positions()[a][i]) + ' '
            b.write(atoms.get_chemical_symbols()[a] + ' ' + pos +
                    ' ' + str(self.hirsh_volrat[a]) + '\n')
        # Write only if periodic
        if atoms.get_pbc().any():
            c = atoms.get_cell()
            for i in range(3):
                b.write('lattice_vector ' + str(c[i][0]) + ' ' +
                        str(c[i][1]) + ' ' + str(c[i][2]) + '\n')
        b.close()

        # Write setting file
        a = open(self.parameters['SETTING_FILE'], 'w+')
        for arg in ['xc',
                    'mbd_cfdm_dip_cutoff',
                    'mbd_scs_dip_cutoff',
                    'mbd_supercell_cutoff',
                    'mbd_scs_vacuum_axis']:
            a.write(arg + ' ' + self.parameters[arg] + '\n')
        a.close()

    def read_results(self):
        """Read energy from output file."""

        a = open(self.parameters['OUTPUT_FILE'], 'r')
        tmp = a.readlines()
        a.close()
        if self.parameters['mode'] == 'TS':
            self.results['energy'] = float(tmp[-2].split()[6])
        elif self.parameters['mode'] == 'rsSCS':
            # In case periodic the output differs
            if self.atoms.get_pbc().any():
                self.results['energy'] = float(tmp[-8].split()[6])
            else:
                self.results['energy'] = float(tmp[-6].split()[6])
        else:
            print('Wohooo no energy for youhoooouu')

    def get_forces(self, apply_constraint=True):
        if self.calculate_forces:
            self.results['forces'] = np.array([[self.calculate_numerical_forces(i, j)
                                                for j in range(3)] for i in range(len(self.atoms))])
            return self.results['forces']
        else:
            return [0.] * self.atoms.get_number_of_atoms()

    def get_potential_energy(self, atoms=None):
        self.calculate(atoms)
        energy = self.get_property('energy', atoms)
        return energy

    def calculate_numerical_forces(self, i, j):
        # based on numeric_force in ase.calculators.test
        delta=0.001
        keep = self.atoms.positions[i, j]
        self.atoms.positions[i, j] += delta
        eplus = self.get_potential_energy(self.atoms)
        self.atoms.positions[i, j] -= 2 * delta
        eminus = self.get_potential_energy(self.atoms)
        self.atoms.positions[i, j] = keep
        return (eminus - eplus) / (2 * delta)

    def set_hirshfeld(self, hirsh_volrat):
        self.hirshvolrat_is_set = True
        self.hirsh_volrat = hirsh_volrat
