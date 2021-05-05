import os 
import torch
import ase
import schnetpack as spk
from ase import Atoms
import numpy as np
import time 
import logging
import argparse 
from sys import argv
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


def get_parser():
    """ Setup parser for command line arguments """
    parser = argparse.ArgumentParser()

    #command-specific
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--cuda', help='Set flag to use GPU(s)', action='store_true')
    parser.add_argument('--nsteps', help='Number of max. optimization steps', type = int, default=400)
    parser.add_argument('--force_mask',help='Atom number of excluded atom for force training', type=int,default=0)
    parser.add_argument('--vdw',help='Define the type of vdw correction. Options: vdw, mbd', default = 'vdw', type = str)
    parser.add_argument('--ts',help='Define the type of vdw correction. Options: TS, TSsurf', default = 'TS', type = str)
    parser.add_argument('--fmax',help='Temperature for Basin Hopping algorithm', default = 0.05, type = float)
    parser.add_argument('--temp',help='Temperature for Basin Hopping algorithm', default = 300, type = int)
    parser.add_argument('initialcondition', help='Input file')
    parser.add_argument('path', help='Path for the simulation')
    parser.add_argument('modelpath', help='Destination for models and logs')
    parser.add_argument('--hirshfeld_modelpath', type=str,help = 'Destination for models and logs for hirshfeld volume rations',default=None)
    parser.add_argument('--overwrite', action='store_true', help='Overwrite old directories')

    return parser

if __name__ == '__main__':

    #parse arguments
    parser = get_parser()
    args = parser.parse_args()

    # Determine the device
    device = torch.device("cuda" if args.cuda else "cpu")

    #read model arguments and load models
    force_model_args =spk.utils.read_from_json(os.path.join(args.modelpath,"args.json"))
    force_model = torch.load(os.path.join(args.modelpath,"best_model"),map_location=device)
    if force_model_args.parallel == True and device == "cpu":
        force_model = force_model.module
    #do the same for the hirshfeld model
    if args.hirshfeld_modelpath is not None:
        hirshfeld_model = torch.load(os.path.join(args.hirshfeld_modelpath,"best_model"),map_location=device)
        hirshfeld_model_args = spk.utils.read_from_json(os.path.join(args.hirshfeld_modelpath,"args.json"))
        if hirshfeld_model_args.parallel == True and device == "cpu": #and args.device == "gpu":
            hirshfeld_model = hirshfeld_model.module
    
    #environment provider for ML
    environment_provider=get_environment_provider(force_model_args,device=device)

    #create or overwrite the folder and change directory
    if args.overwrite and os.path.exists(args.path):
        logging.info('Existing folder will be overwritten...')
        rmtree(args.path)
    if not os.path.exists(args.path):
        os.makedirs(args.path)
    current_path = os.getcwd()
    os.chdir(args.path)

    #read atoms object
    atomspath = os.path.join(current_path,args.initialcondition)
    atoms_init = ase.io.read(os.path.join(current_path,args.initialcondition))
    
    natoms = atoms_init.get_number_of_atoms()
    if atoms_init.get_pbc()[0]==True:
        stress = spk.Properties.stress
        force_model.requires_stress = True

    else:
        stress = None
    qm_calc = spk_vdw_interface.SpkVdwCalculator(force_model,
                                                 hirshfeld_model = hirshfeld_model,
                                                 device = device,
                                                 energy = spk.Properties.energy,
                                                 forces = spk.Properties.forces,
                                                 hirsh_volrat = "hirshfeld_volumes",
                                                 energy_units = 'eV', forces_units='eV/A',
                                                 environment_provider = environment_provider)



    vdw_calc = MBD(
        scheme=args.vdw, #VDW or MBD
        params=args.ts,  #TS or TSsurf
        ts_sr=0.94,  #for vdw
        #beta = 0.83 #for MBD
        k_grid=(4,4,1))
    
    dispcorr = DispersionCorrectionCalculator(
                qm_calculator = qm_calc,
                mm_calculator = vdw_calc,
                )

    atoms_init.set_calculator(dispcorr)
    e = atoms_init.get_potential_energy()
    f = atoms_init.get_forces()
    opt = BFGS(atoms_init,trajectory='opt.traj')
    opt.run(fmax=args.fmax)
