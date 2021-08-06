This is a repository of a SchNetPack Calculator capable to predict energies, forces, and hirshfeld volume rations that can subsequently be used to correct energies and forces using many-body dispersion correction or van-der-Waals correction.
It is based on the Atomic Simulation Environment and SchNetPack. 

Requirements:
1. SchNetPack
2. ase
3. libmbd

The following installation is recommended:

1. download anaconda (https://docs.anaconda.com/)
2. create a new environmet: conda create -n ml python
3. The environment is named "ml", activate it via: conda activate ml
4. To install libmbd, execute the following:
    a. conda install pip
    b. pip install --upgrade pip setuptools wheel
    c. conda install -c conda-forge libmbd
    d. conda install gcc_linux-64
    e. conda install gxx_linux-64
    f. pip install pymbd
5. Install SchNetPack including all requirements (for the correct pytorch cuda combination, check: https://pytorch.org/get-started/locally/). A suggested installation procedure could be:
    a. conda install ase numpy six scipy matplotlib  h5py tqdm pytest pytest-datadir
    b. conda install pytorch torchvision torchaudio cudatoolkit=XX.X -c pytorch
    c. see https://github.com/atomistic-machine-learning/schnetpack for installation of SPK
6. Install this github repository: execute "pip install ." in this directory


Once done, don't forget to set your LD_LIBRARY_PATH.

.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _IRC: http://webchat.freenode.net/?randomnick=0&channels=ase
