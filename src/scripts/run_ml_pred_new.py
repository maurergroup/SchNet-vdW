import os 
import torch
import ase
import argparse
import schnetpack as spk
from ase import Atoms
from sys import argv
import numpy as np
import time 
import logging
from shutil import copyfile, rmtree
from ase import units
from ase.optimize import QuasiNewton
from schnetpack.utils.script_utils.settings import get_environment_provider
from schnetpack.utils.script_utils.parsing import build_parser
from spk_vdw.mbd import MBD
from spk_vdw import qmme
from spk_vdw import spk_vdw_interface
from ase.io import read
from ase.optimize import BFGS
from ase.units import kB
from ase.optimize.basin import BasinHopping
from ase.constraints import FixAtoms


def get_parser():
    """ Setup parser for command line arguments """
    parser = argparse.ArgumentParser()

    #command-specific
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--cuda', help='Set flag to use GPU(s)', action='store_true')
    parser.add_argument('--nsteps', help='Number of max. optimization steps', type = int, default=400)
    parser.add_argument('--force_mask',help='Atom number of excluded atom for force training', type=int,default=0)
    parser.add_argument('--mode', help='Type of calculation. Options: sp, opt, md, opt_vdw, bh',default='opt_vdw',type=str)
    parser.add_argument('--vdw',help='Define the type of vdw correction. Options: vdw, mbd', default = 'vdw', type = str)
    parser.add_argument('--fmax',help='Temperature for Basin Hopping algorithm', default = 0.05, type = float)
    parser.add_argument('--temp',help='Temperature for Basin Hopping algorithm', default = 300, type = int)
    parser.add_argument('initialcondition', help='Input file')
    parser.add_argument('path', help='Path for the simulation')
    parser.add_argument('modelpath', help='Destination for models and logs')
    parser.add_argument('--hirshfeld_modelpath', type=str,help = 'Destination for models and logs for hirshfeld volume rations',default=None)
    parser.add_argument('--overwrite', action='store_true', help='Overwrite old directories')

    return parser
def get_initial_position(path):
    
    #blabla
    initfile = open(path,"r").readlines()
    cell=np.zeros((3,3))
    #for i in range(3):
    #    for j in range(3):
    #        cell[i][j]=float(initfile[i].split()[j])
    all_atoms=[]
    atypes=[]
    for i in range(2,len(initfile)):
        atoms=[]
        atypes.append(str(initfile[i].split()[0]))
        for xyz in range(3):
            atoms.append(float(initfile[i].split()[xyz+1]))
        all_atoms.append(np.array(atoms))
    all_atoms=np.array(all_atoms)
    atypes=np.array(atypes)
    atom_object = Atoms(atypes,all_atoms)
    #atom_object.set_cell(cell)
    #atom_object.set_pbc(True)
    #print(atom_object)
    #print(atom_object.get_positions().shape)

    return atom_object


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
 
    #create or overwrite the folder and change directory
    if args.overwrite and os.path.exists(args.path):
        logging.info('Existing folder will be overwritten...')
        rmtree(args.path)
    if not os.path.exists(args.path):
        os.makedirs(args.path)
    current_path = os.getcwd()
    os.chdir(args.path)

    #read initial atom file
    atomspath = os.path.join(current_path,args.initialcondition)
    atoms_init = ase.io.read(os.path.join(current_path,args.initialcondition))
    
    #get environment provider for SchNet models
    environment_provider=get_environment_provider(force_model_args,device=device)
   
    # do a single point calculations
    if args.mode == "sp":
        calculator = spk.interfaces.SpkCalculator(model = force_model, 
                                                  device=device, energy = spk.Properties.energy, 
                                                  forces=spk.Properties.forces,
                                                  force_mask = args.force_mask,
                                                  energy_units='eV', 
                                                  forces_units="eV/A",environment_provider=environment_provider)
        atoms_init.set_calculator(calculator)
        print("Prediction:", atoms_init.get_total_energy(), atoms_init.get_forces())
    
    # do an optimization of MD run
    if args.mode == "opt" or args.mode =="md":
        #set a force mask to mask all values that should stay constant
        model_ase = spk_vdw_interface.SpkVdwCalculator(atomspath,force_model,
                                                os.getcwd(),device=device,
                                                energy=spk.Properties.energy,
                                                forces=spk.Properties.forces,
                                                energy_units='eV', forces_units="eV/A",
                                                force_mask=args.force_mask,environment_provider=environment_provider)
        if args.mode == "opt":
            model_ase.optimize(fmax=args.fmax,steps=args.nsteps)    
            model_ase.compute_normal_modes()

        if args.mode == "md":
            T=open("md_input","r").readlines()
            for line in T:
                if line.startswith("Temp"):
                    temp=int(line.split()[1])
                if line.startswith("time"):
                    time=int(line.split()[1])
            model_ase.init_md('simulation',temp_bath=temp,reset=True)
            model_ase.run_md(time)
            results=np.loadtxt(os.path.join(os.getcwd(),'simulation.log'),skiprows=1)

    # make an optimization with vdw correction 
    if args.mode =="opt_vdw" or args.mode =="bh":
        
        natoms = atoms_init.get_number_of_atoms()
        #atoms_init.set_cell(([1,0,0],[0,1,0],[0,0,1]))
        #print(atoms_init.cell)
        #calculator = spk.interfaces.AseVdwInterface(atomspath,force_model, os.getcwd(),device=device,
        #                                         hirshfeld_model = hirshfeld_model, 
        #                                         energy = spk.Properties.energy, 
        #                                         forces=spk.Properties.forces,
        #                                         force_mask = args.force_mask,
        #                                         hirshfeld="hirshfeld_volumes",
        #                                         energy_units='eV', forces_units="eV/A",
        #                                         environment_provider=environment_provider)
        calculator = spk_vdw_interface.SpkVdwCalculator(force_model, 
                                                     hirshfeld_model = hirshfeld_model,
                                                     device = device,
                                                     energy = spk.Properties.energy,
                                                     forces = spk.Properties.forces,
                                                     # set in qmme force_mask = args.force_mask,
                                                     hirsh_volrat = "hirshfeld_volumes",  
                                                     energy_units = 'eV', forces_units='eV/A',
                                                     environment_provider = environment_provider)
        if args.vdw == "vdw":
            vdw_calc = MBD(
                scheme="VDW", 
                params="TS", #TODO currently hardcoded, should be in args
                ts_sr=0.94, #TODO currently hardcoded to PBE value, should be in args
                k_grid=None) #TODO currently hardcoded, but needs to be passed in depending on aperiodic or PBC

        if args.vdw == "mbd":
            vdw_calc = MBD(
                scheme="MBD", 
                params="TSSURF", #TODO currently hardcoded, should be in args
                beta=0.83, #TODO currently hardcoded to PBE value, should be in args
                k_grid=None)
        #atoms_init.set_constraint(FixAtoms(indices=np.arange(args.force_mask)))
        #print(atoms_init)
        VDW = qmme.qmme(atoms=atoms_init,
                nqm_regions = 1,
                nmm_regions = 1,
                qm_pbcs = [True,True,False],
                mm_pbcs = [True,True,False],
                qm_calculators = [calculator],
                mm_calculators = [vdw_calc],
                qm_atoms = [[(0,natoms)]],
                mm_cell = [np.zeros((3,3))],
                qm_cell = [np.zeros((3,3))],
                mm_atoms = [[(0,natoms)]],
                freeze = np.arange(args.force_mask),
                mm_mode = "explicit")

        atoms_init.set_calculator(VDW)
        #model_ase.optimize(fmax=args.fmax,steps=100)
        filename  = "opt"
        
        if args.mode == "opt_vdw":
    
            optimizer = BFGS( atoms_init, trajectory = "%s.traj" %filename, restart="%s.pkl" %filename)
            optimizer.run(fmax=args.fmax,steps=args.nsteps)
    
        if args.mode == "bh":

            # Set up the basin hopping system
            bh = BasinHopping(atoms=atoms_init,         # the system to optimize
                      temperature=300 * kB, # 'temperature' to overcome barriers
                      dr=0.5,               # maximal stepwidth
                      optimizer=BFGS,      # optimizer to find local minima
                      fmax=args.fmax,             # maximal force for the optimizer
                      )

            bh.run(100)
