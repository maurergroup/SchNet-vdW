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
from spk_vdw import spk_vdw_interface, spk_vdw_qbc_interface
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
    parser.add_argument('--mode',help='Define the type of calculations. Options: optimization (key: opt), basin hopping (key: bh)', default = 'opt', type = str)
    parser.add_argument('--dr', help = "Maximal stepwidth in basin hopping.", default = 0.5, type = float)
    parser.add_argument('--adaptive_fmax', help = "Fmax at which energy variance is included in optimization (early stopping).", default = 0.15, type = float)
    parser.add_argument('initialcondition', help='Input file')
    parser.add_argument('--qbc', help="Enable Querry by committee", action='store_true')
    parser.add_argument('path', help='Path for the simulation')
    parser.add_argument('modelpath', help='Destination for models and logs')
    parser.add_argument('--Hessianfile', help = 'Specify the filename that contains the Hessian to initialize geometry optimizations', type = str)
    parser.add_argument('--hirshfeld_modelpath', type=str,help = 'Destination for models and logs for hirshfeld volume rations',default=None)
    parser.add_argument('--overwrite', action='store_true', help='Overwrite old directories')
    parser.add_argument('--extrapolate', action='store_true', help='Extrapolative system not included in the training set (size, molecular system).')
    parser.add_argument('--kgrid', type = list, nargs = 3, help = "Specify k-grid (x,y,z)")
    parser.add_argument('--nmodels', help = "specify the number of ML models to use", default = 1, type = int)
    parser.add_argument('--nhmodels', help = "specify the number of Hirshfeld models to use", default = 1, type = int)
    parser.add_argument('--bhsteps', help = "Number of basin hopping steps", default = 100, type = int)
    return parser

if __name__ == '__main__':

    #parse arguments
    parser = get_parser()
    args = parser.parse_args()

    # Determine the device
    device ="cpu"# torch.device("cuda" if args.cuda else "cpu")

    #read model arguments and load models
    #modelpath="MLModels/Au_C/EF/"
    #hmodel = "MLModels/Au_C/Hirshfeld/"
    force_model_args={}
    kgrid =[]
    if args.kgrid is not None:
        kgrid.append(int(args.kgrid[0][0]))
        kgrid.append(int(args.kgrid[1][0]))
        kgrid.append(int(args.kgrid[2][0]))
    else:
        kgrid = None
    force_model = {}
    hirshfeld_model = {}
    hirshfeld_model_args = {}
    environment_provider = {}
    for i in range(args.nmodels):
        force_model_args[i] = spk.utils.read_from_json(os.path.join(args.modelpath+"/Model%i/"%(i+1),"args.json"))
        force_model[i] = torch.load(os.path.join(args.modelpath+"/Model%i/"%(i+1),"best_model"),map_location=device)
        if force_model_args[i].parallel == True and device == "cpu":
            force_model[i] = force_model[i].module
        #environment provider for ML, 
        environment_provider[i]=get_environment_provider(force_model_args[i],device=device)

    #do the same for the hirshfeld model
    for i in range(args.nhmodels):
        hirshfeld_model[i] = torch.load(os.path.join(args.hirshfeld_modelpath+"/Model%i/"%(i+1),"best_model"),map_location=device)
        hirshfeld_model_args[i] = spk.utils.read_from_json(os.path.join(args.hirshfeld_modelpath+"/Model%i/"%(i+1),"args.json"))
        if hirshfeld_model_args[i].parallel == True and device == "cpu": #and args.device == "gpu":
            hirshfeld_model[i] = hirshfeld_model[i].module
    
    #read atoms object
    atoms_init = ase.io.read(args.initialcondition)
    import ase.constraints
    
    natoms = atoms_init.get_number_of_atoms()
    if atoms_init.get_pbc()[0]==True:
        stress = spk.Properties.stress
        for i in range(args.nmodels):
            force_model[i].requires_stress = True

    else:
        stress = None
    qm_calc = spk_vdw_qbc_interface.SpkVdwCalculator(force_model,
                                                 hirshfeld_model = hirshfeld_model,
                                                 device = device,
                                                 energy = spk.Properties.energy,
                                                 forces = spk.Properties.forces,
                                                 hirsh_volrat = "hirshfeld_volumes",
                                                 energy_units = 'eV', forces_units='eV/A',
                                                 environment_provider = environment_provider,
                                                 fmax=args.fmax,
                                                 nmodels = args.nmodels,
                                                 nhmodels = args.nhmodels,
                                                 extrapolate = args.extrapolate,
                                                 adaptive_fmax = args.adaptive_fmax,
                                                 qbc = args.qbc,
                                                 logfile = args.path + "/opt.log",
                                                 mode=args.mode)
 
    vdw_calc = MBD(
        scheme=args.vdw, #VDW or MBD
        params=args.ts,  #TS or TSsurf
        ts_sr=0.94,  #for vdw
        #beta = 0.83 #for MBD
        k_grid=kgrid)
    
    dispcorr = DispersionCorrectionCalculator(
                qm_calculator = qm_calc,
                mm_calculator = vdw_calc,
                )

    atoms_init.set_calculator(dispcorr)
    #c=ase.constraints.FixAtoms(indices=[atom.index for atom in atoms_init if atom.symbol == 'C'] )
    #atoms_init.set_constraint(c) #ase.constraints.FixAtoms(np.arange(686)))
    if args.mode == "opt":
        opt = BFGS(atoms_init,trajectory=args.path+'/opt.traj')
        if args.Hessianfile == True:
            opt.H0 = np.loadtxt(args.Hessianfile)
        opt.run(fmax=args.fmax)
    elif args.mode == "bh":
        bh = BasinHopping( atoms=atoms_init,
             temperature = args.temp * kB, # Temperature to overcome barriers
             dr = args.dr, # maximal stepwidth
             optimizer = BFGS, #optimizer to find local minima
             fmax = args.fmax, #maximal force for the optimizer
             trajectory = args.path + '/lowest.traj',
             local_minima_trajectory = args.path + '/local_minima.traj'
             )
        bh.run(args.bhsteps)
    else:
        print("Mode not implemented.") 
