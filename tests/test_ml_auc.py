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
    #parser = get_parser()
    #args = parser.parse_args()

    # Determine the device
    device ="cpu"# torch.device("cuda" if args.cuda else "cpu")

    #read model arguments and load models
    modelpath="MLModels/Au_C/EF/"
    hmodel = "MLModels/Au_C/H/"
    force_model_args =spk.utils.read_from_json(os.path.join(modelpath,"args.json"))
    force_model = torch.load(os.path.join(modelpath,"best_model"),map_location=device)
    if force_model_args.parallel == True and device == "cpu":
        force_model = force_model.module
    #do the same for the hirshfeld model
    hirshfeld_model = torch.load(os.path.join(hmodel,"best_model"),map_location=device)
    hirshfeld_model_args = spk.utils.read_from_json(os.path.join(hmodel,"args.json"))
    if hirshfeld_model_args.parallel == True and device == "cpu": #and args.device == "gpu":
        hirshfeld_model = hirshfeld_model.module
    
    #environment provider for ML
    environment_provider=get_environment_provider(force_model_args,device=device)

    #read atoms object
    atoms_init = ase.io.read("50.in")
    import ase.constraints
    
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
						 stress = stress,
                                                 hirsh_volrat = "hirshfeld_volumes",
                                                 energy_units = 'eV', forces_units='eV/A',
                                                 environment_provider = environment_provider)



    vdw_calc = MBD(
        scheme="VDW", #VDW or MBD
        params="TS",  #TS or TSsurf
        ts_sr=0.94,  #for vdw
        #beta = 0.83 #for MBD
        k_grid=None)#(4,4,1))
    
    dispcorr = DispersionCorrectionCalculator(
                qm_calculator = qm_calc,
                mm_calculator = vdw_calc,
                )

    atoms_init.set_calculator(dispcorr)
    c=ase.constraints.FixAtoms(indices=[atom.index for atom in atoms_init if atom.symbol == 'C'] )
    atoms_init.set_constraint(c) #ase.constraints.FixAtoms(np.arange(686)))
    e = atoms_init.get_potential_energy()
    f = atoms_init.get_forces()
    opt = BFGS(atoms_init,trajectory='opt.traj')
    opt.run(fmax=0.05)
