import os 
import torch
import ase
import schnetpack as spk
from ase import Atoms
import numpy as np
import time 
import logging
from shutil import copyfile, rmtree
from ase import units
from ase.optimize import BFGS
from schnetpack.utils.script_utils.settings import get_environment_provider
from schnetpack.utils.script_utils.parsing import build_parser
from spk_vdw.mbd import MBD
from spk_vdw.dispcorr import DispersionCorrectionCalculator
from spk_vdw import spk_vdw_interface
from ase.io import read
from ase.units import kB
from ase.optimize.basin import BasinHopping
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.calculators.aims import Aims

ang = 1.0

def benzene_dimer():
    return [
        (np.array([
            (-1.047, -1.421, 0.000), (-1.454, -0.855, 1.206),
            (-1.454, -0.855, -1.206), (-2.266, 0.277, 1.206),
            (-2.671, 0.845, 0.000), (-2.266, 0.277, -1.206),
            (-1.133, -1.292, -2.142), (-2.582, 0.716, -2.143),
            (-3.303, 1.723, 0.000), (-2.582, 0.716, 2.143),
            (-1.133, -1.292, 2.142), (-0.406, -2.291, 0.000)]) * ang,
            6 * ['C'] + 6 * ['H'],
            [0.825, 0.821, 0.821, 0.815, 0.814, 0.815,
             0.624, 0.611, 0.610, 0.611, 0.624, 0.643]),
        (np.array([
            (1.047, 1.421, 0.000), (1.454, 0.855, -1.206),
            (1.454, 0.855, 1.206), (2.266, -0.277, -1.206),
            (2.671, -0.845, 0.000), (2.266, -0.277, 1.206),
            (0.406, 2.291, 0.000), (1.133, 1.292, 2.142),
            (2.582, -0.716, 2.143), (3.303, -1.723, 0.000),
            (2.582, -0.716, -2.143), (1.133, 1.292, -2.142)]) * ang,
            6 * ['C'] + 6 * ['H'],
            [0.825, 0.821, 0.821, 0.815, 0.814, 0.815,
             0.643, 0.624, 0.611, 0.610, 0.611, 0.624])]

mol1, mol2 = benzene_dimer()
benzene1 = Atoms(mol1[1], positions=mol1[0])
benzene2 = Atoms(mol2[1], positions=mol2[0])

dimer = benzene1 + benzene2
natoms = len(dimer)

qm_calc = Aims(
    aims_command='/home/chem/mssdjc/software/fhiaims/code/binaries/SCRTP/aims.180128.mpi.x',
    outfilename = "aims.out",
    species_dir='/home/chem/mssdjc/software/fhiaims/species_defaults/light',
    xc='PBE',
    occupation_type = ['gaussian',0.01],
    sc_iter_limit = 100,
    #spin = 'collinear',
    relativistic = ['atomic_zora','scalar'],
    sc_accuracy_etot=1e-6,
    sc_accuracy_eev=1e-3,
    sc_accuracy_rho=1e-6,
    sc_accuracy_forces=1e-4,
    #load_balancing = True,
    #k_grid = None,
    #restart_aims='wvfn.dat',
    output = ['hirshfeld'],
)

vdw_calc = MBD(
    scheme="VDW", #VDW or MBD
    params="TS",  #TS or TSSURF
    ts_sr=0.94,  #for vdw
    #beta = 0.83 #for MBD
    k_grid=None)

dispcorr = DispersionCorrectionCalculator(
                qm_calculator = qm_calc,
                mm_calculator = vdw_calc,
                )

dimer.set_calculator(dispcorr)

e = dimer.get_potential_energy()
f = dimer.get_forces()
print(e)
print(f)

opt = BFGS(dimer,trajectory='opt.traj')
opt.run(fmax=0.05)
