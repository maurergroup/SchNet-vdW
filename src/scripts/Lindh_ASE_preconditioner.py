#This script was provided by Dimitrii Maksimov
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.lj import LennardJones
import sys
import os


atoms = read('geometry.in', format = 'aims')
calculator = LennardJones()   
atoms.set_calculator(calculator)
opt = BFGS(atoms)

#------------------------------------------------------
# Create Lindh Hessian matrix
os.system('python2.7 Lindh.py geometry.in --full Hessian.dat') 
# Load Initial Hessian
opt.H0 = np.loadtxt('Hessian.dat') 
#------------------------------------------------------

opt.run(fmax=1e-3, steps=300)

