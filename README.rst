This is a repository of a SchNetPack Calculator capable to predict energies, forces, and Hirshfeld volume ratios that can subsequently be used to correct energies and forces using many-body dispersion correction or van der Waals correction. It is based on the Atomic Simulation Environment (ASE) and SchNetPack. 

Requirements:
1. SchNetPack
2. ASE
3. Libmbd

The following installation procedure is recommended:

1. Download Anaconda (https://docs.anaconda.com/)
2. Create a new environment: conda create -n ml python
3. The environment is named "ml", activate it via: conda activate ml
4. To install Libmbd, execute the following:
    a. conda install pip
    b. pip install --upgrade pip setuptools wheel
    c. conda install -c conda-forge libmbd
    d. conda install gcc_linux-64
    e. conda install gxx_linux-64
    f. pip install pymbd
5. Install SchNetPack including all requirements (for the correct PyTorch CUDA combination, see: https://pytorch.org/get-started/locally/). A suggested installation procedure is:
    a. Clone SchNetPack (https://github.com/atomistic-machine-learning/schnetpack) 
    b. Go into the SchNetPack folder and install it via "pip install ."
Note that installation does not work via Conda and will report conflicts.
6. Install this (https://github.com/maurergroup/SchNet-vdW/) github repository: execute "pip install ." in this directory

Once done, don't forget to set your LD_LIBRARY_PATH.

.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _IRC: http://webchat.freenode.net/?randomnick=0&channels=ase
