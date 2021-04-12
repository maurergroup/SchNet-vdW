#!/bin/bash
#SBATCH -J BH20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=48:00:00
ulimit -s unlimited

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

MY_NUM_THREADS=$SLURM_CPUS_PER_TASK

CWD=/home/chem/mssdjc/Research/GlobalOpt/TestSetOpt/P2O/1
#do optimizations
MODEL1=/home/chem/mssdjc/Research/GlobalOpt/Training/EF
MODEL2=/home/chem/mssdjc/Research/GlobalOpt/Training/Hirshfeld
cd $CWD
source /home/chem/mssdjc/software/anaconda3/etc/profile.d/conda.sh

python $VDW/run_ml_pred.py geometry.in Opt $MODEL1 --hirshfeld_modelpath $MODEL2 --mode opt_vdw --vdw dftdisp --fmax 0.05  >> opt.log

