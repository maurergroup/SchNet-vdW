This is a repository of a SchNetPack Calculator capable to predict energies, forces, and Hirshfeld volume ratios that can subsequently be used to correct energies and forces using many-body dispersion correction or van der Waals correction. It is based on the Atomic Simulation Environment (ASE) and SchNetPack. 

Requirements:
1. SchNetPack
2. ASE
3. Libmbd

The following installation procedure is recommended when using anaconda. See below for an installation on HPE Cray Supercomputers:

1. Download Anaconda (https://docs.anaconda.com/) and install it.
2. Make sure you can create environments: ``conda config --append channels conda-forge``
4. To install Libmbd, several packages have to be installed before, otherwise an error might be reported. For us, the following procedure worked:
    a.  ``conda create -n vdw python h5py tensorboardX pytorch ase numpy six protobuf scipy matplotlib python-dateutil pyyaml tqdm pyparsing kiwisolver cycler netcdf4 hdf5 h5utils jupyter -c openbabel``
    b. The environment is named "vdw", activate it via: ``conda activate vdw``
    c. ``conda install pip``
    d. ``conda install -c conda-forge libmbd``
    e. ``pip install pymbd``
5. Install SchNetPack including all requirements (for the correct PyTorch CUDA combination, see: https://pytorch.org/get-started/locally/). A suggested installation procedure is:
    a. Clone SchNetPack (https://github.com/juliawestermayr/schnetpack.git or https://github.com/atomistic-machine-learning/schnetpack) 
    b. Go into the SchNetPack folder and install it via ``python setup.py install``
Note that installation does not work via Conda and will report conflicts.
6. Install this (https://github.com/maurergroup/SchNet-vdW/) github repository: execute ``python setup.py install `` in this directory

Once done, don't forget to set your LD_LIBRARY_PATH.

.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _IRC: http://webchat.freenode.net/?randomnick=0&channels=ase


Installing Libmbd on Cray:

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
