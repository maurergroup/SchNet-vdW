This is a repository of a SchNetPack Calculator capable to predict energies, forces, and Hirshfeld volume ratios that can subsequently be used to correct energies and forces using many-body dispersion correction or van der Waals correction. It is based on the Atomic Simulation Environment (ASE) and SchNetPack. 
Installation with anaconda is recommended.

Requirements:
1. SchNetPack
2. ASE
3. Libmbd

The following installation procedure is recommended when using anaconda:

1. Download Anaconda (https://docs.anaconda.com/) and install the latest version (we used 3.9.7 when creating this README file):
    i. ``wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh``
    ii. ``sh Anaconda3-2021.11-Linux-x86_64.sh``
    
2. Make sure you can create environments: ``conda config --append channels conda-forge``

3. Install pymbd (https://pypi.org/project/pymbd/)
    i. ``conda install -c conda-forge libmbd``
    ii. Install pip if you don't have it already:  ``conda install pip``
    iii. ``pip install pymbd``

4. Install SchNetPack including all requirements (for the correct PyTorch CUDA combination, see: https://pytorch.org/get-started/locally/). A suggested installation procedure is: ``pip install schnetpack``

5. ML for hirshfeld volume ratios:
    i. Download the code to use ML models for Hirshfeld volume ratios: ``wget https://figshare.com/ndownloader/files/34002485?private_link=78b54de875cfb9cadbdd``
    ii. ``unzip 34002485?private_link=78b54de875cfb9cadbdd``
    iii. Go into the SchNet_EV-AuC_stable folder.
    iv. Install it via ``python setup.py install``

6. Install SchNet-vdW:
    i. Go into the SchNet-vdW folder.
    ii. Install it via ``python setup.py install``
    iii. Go into the ``src/spk_vdw`` folder and compile the code by executing the commands written in the files "build_options_sdc_gnu" or "build_options_sdc".
    iv. This should generate a file "sdc.cpython-...-gnu.so". Copy this file to "sdc.so".

7. Set your LD_LIBRARY_PATH to the lib-folder of your anaconda environment. e.g. ``export LD_LIBRARY_PATH=/home/software/anaconda3/envs/vdw/lib``


Note that installation does not work via Conda and will report conflicts.

.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _IRC: http://webchat.freenode.net/?randomnick=0&channels=ase


GETTING STARTED:

1. Go into the SchNet-vdW/tests folder.

2a. Running DFT optimizations: Note that you need fhi-aims for this purpose
   Execute one of the available python scripts that contain "mbd" in their names.
   
2b. Running ML optimizations:

   Some example commands can be found in the "test_ml.sh" file.
   
   ML models are provided in the "MLModels" folder. Execute the following command for a quick test:
   
   ``python test_ml.py B2O.in OptX2O MLModels/X2O/EF --hirshfeld_modelpath MLModels/X2O/Hirshfeld/ --vdw vdw --ts TSSurf --Hessianfile OptX2O/HessianB2O.dat --kgrid 7 7 1 >      OptX2O/opt.log``
   
   The input file is called B2O.in. The files will be saved in the OptX2O-folder. The path of the ML models for energies and forces. The path for hirshfeld models can be        indicated via --hirshfeld_modelpath; --vdw specifies the type of vdW-correction. --ts can be used to specify TS or TSSurf. The Hessianfile-flag can be set to specify an initial Hessian guess. --kgrid can be used in case a k-grid is needed for calculations. The output will be saved in the "opt.log" file in the corresponding folder.
   
   Hessian input: In case you want to initialize optimizations with a Lindh-Hessian, you can use the Lindh.py file provided by within FHI-aims (https://gitlab.com/users/sign_in). Alternatively, you can use Lindh-Hessian generation scripts from the gensec repository (https://github.com/sabia-group/gensec).
    



Example for installing Libmbd on Cray:

1. add the following lines to your ~/.profile:

``module load PrgEnv-gnu``

``module load cray-python``

``module load cmake``

``export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/opt/cray/xpmem/2.2.40-7.0.1.0_2.7__g1d7a24d.shasta/lib64/pkgconfig``


2. install libmbd via:
    a. ``wget https://github.com/libmbd/libmbd/releases/download/0.12.3/libmbd-0.12.3.tar.gz``
    b. ``tar -xf libmbd-0.12.3.tar.gz``
    c. ``python -m venv venv --system-site-packages``
    d. ``source venv/bin/activate``
    e. ``cmake -S libmbd-0.12.3 -B build -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_Fortran_FLAGS="-fPIC -O3 -fallow-argument-mismatch" -DENABLE_SCALAPACK_MPI=ON -DSCALAPACK_LIBRARY= -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV``
    f. ``make -C build install``
    g. ``env LIBMBD_PREFIX=$VIRTUAL_ENV pip install pymbd[mpi,test]==0.12.3``

3. you can test it by running: 
``time srun --cpu_bind=cores --distribution=block:cyclic --hint=nomultithread python -m pymbd.benchmark --finite``

4. Upgrade pip via:
``pip install --upgrade pip setuptools wheel``

5. Install SchNetPack and Spk-vdw (this repository) as outlined above.
