import os
import io

from setuptools import setup, find_packages


def read(fname):
    with io.open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8') as f:
        return f.read()

setup(
    name='spk_vdw',
    version='0.0.0',
    author='',
    email='julia.westermayr@warwick.ac.uk',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    scripts=['./scripts/run_ml_pred_new.py','./scripts/run_ml_pred.py','./spk_vdw/dftdisp.py','./spk_vdw/dispcorr.py','./spk_vdw/mbd.py','./spk_vdw/spk_vdw_qbc_interface.py','./spk_vdw/spk_vdw_interface.py','./spk_vdw/__init__.py','./spk_vdw/mbdio.py','./spk_vdw/qmme.py']
    python_requires='>=3.6',
    install_requires=[
        "torch>=0.4.1",
        "numpy",
        "ase>=3.16",
        "tensorboardX",
        "h5py"
    ]
)
