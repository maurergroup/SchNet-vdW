This is a repository of a SchNetPack calculator capable of predicting energies, forces, and Hirshfeld volume ratios that can subsequently be used to correct energies and forces using a many-body dispersion correction or a van der Waals correction. It is based on the Atomic Simulation Environment (ASE) and SchNetPack. 
Installation with Anaconda is recommended.

Installation Procedure
=======================

Requirements:
1. SchNetPack
2. ASE
3. Libmbd (attention: only works on Mac or Linux)

1. Download Anaconda (https://docs.anaconda.com/) and install the latest version (we used 3.9.7 when creating this README file): ``wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh`` and then ``sh Anaconda3-2021.11-Linux-x86_64.sh`` 

2. Make sure you can create environments: ``conda config --append channels conda-forge``

3. Install pymbd (https://pypi.org/project/pymbd/): ``conda install -c conda-forge libmbd``. Install pip if you don't have it (``conda install pip``) and finally ``pip install pymbd`` Attention: please install pymbd before schnetpack. If it's done the other way around, conflicts may arise.

4. Install SchNetPack using this repository: https://github.com/juliawestermayr/schnetpack (v0.1) and install all requirements (``pip install -r requirements.txt`` in the ``schnetpack`` directory) to train on hirshfeld volume ratios. A tutorial is provided on figshare: 10.6084/m9.figshare.19134602
If you don't want to follow the tutorial, but just want to do some optimizations with the existing code and ML models, please download the code to use ML models for Hirshfeld volume ratios: ``wget https://figshare.com/ndownloader/files/34002485?private_link=78b54de875cfb9cadbdd``. After this: ``unzip 34002485?private_link=78b54de875cfb9cadbdd``. After changing into the ``SchNet_EV-AuC_stable`` directory, install it via ``python setup.py install``

6. Install SchNet-vdW by changing into the ``SchNet-vdW`` directory and running ``python setup.py install``

7. Set your ``LD_LIBRARY_PATH`` variable to the ``lib`` directory of your Anaconda environment. e.g. ``export LD_LIBRARY_PATH=/home/software/anaconda3/envs/vdw/lib``


.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _IRC: http://webchat.freenode.net/?randomnick=0&channels=ase


Getting Started:
================

We provide a short tutorial including training of ML models on figshare (10.6084/m9.figshare.19134602), if you are not interested and want to test optimization, please follow the procedure below.

1. Change into the ``SchNet-vdW/tests`` directory.

2a. Running DFT optimizations: Note that you need FHI-aims for this purpose. Execute one of the available Python scripts containing "mbd" in their names.

2b. Running ML optimizations: Some example commands can be found in the ``test_ml.sh`` file. ML models are provided in the ``MLModels`` folder. Execute the following command for a quick test:
   
   ``python test_ml.py B2O.in OptX2O MLModels/X2O/EF --hirshfeld_modelpath MLModels/X2O/Hirshfeld/ --vdw vdw --ts TSSurf --Hessianfile OptX2O/HessianB2O.dat --kgrid 7 7 1 >      OptX2O/opt.log``
   
   The input file is called B2O.in. The files will be saved in the OptX2O-folder. The path of the ML models for energies and forces. The path for Hirshfeld models can be indicated via ``--hirshfeld_modelpath``; ``--vdw`` specifies the type of vdW-correction. ``--ts`` can be used to specify TS or TSSurf. The ``--Hessianfile`` flag can be set to specify an initial Hessian guess. ``--kgrid`` can be used in case a k-grid is needed for calculations. The output will be saved in the ``opt.log`` file in the corresponding folder.
   
   Hessian input: In case you want to initialize optimizations with a Lindh-Hessian, you can use the ``Lindh.py`` file provided by within FHI-aims (https://gitlab.com/users/sign_in). Alternatively, you can use Lindh-Hessian generation scripts from the gensec repository (https://github.com/sabia-group/gensec).
    
Example for installing Libmbd on Cray:
=======================================

1. Add the following lines to your ``~/.profile`` file:
    a. ``module load PrgEnv-gnu``
    b. ``module load cray-python``
    c. ``module load cmake``
    d. ``export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/opt/cray/xpmem/2.2.40-7.0.1.0_2.7__g1d7a24d.shasta/lib64/pkgconfig``
    
2. Install Libmbd via
    a. ``wget https://github.com/libmbd/libmbd/releases/download/0.12.3/libmbd-0.12.3.tar.gz``
    b. ``tar -xf libmbd-0.12.3.tar.gz``
    c. ``python -m venv venv --system-site-packages``
    d. ``source venv/bin/activate``
    e. ``cmake -S libmbd-0.12.3 -B build -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_Fortran_FLAGS="-fPIC -O3 -fallow-argument-mismatch" -DENABLE_SCALAPACK_MPI=ON -DSCALAPACK_LIBRARY= -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV``
    f. ``make -C build install``
    g. ``env LIBMBD_PREFIX=$VIRTUAL_ENV pip install pymbd[mpi,test]==0.12.3``

3. The installation can be tested by running: ``time srun --cpu_bind=cores --distribution=block:cyclic --hint=nomultithread python -m pymbd.benchmark --finite``

4. Upgrade pip via ``pip install --upgrade pip setuptools wheel``

5. Install SchNetPack and Spk-vdw (this repository) as outlined above.


The tagged version v0.1 was used for arXiv:2202.13009v1 and for tutorials provided on figshare: https://figshare.com/articles/software/SchNet_vdW/19134602
