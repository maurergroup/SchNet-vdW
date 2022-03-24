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
    scripts=[],
    python_requires='>=3.6',
    install_requires=[
        "torch>=1",
        "numpy",
        "ase>=3.16",
        "h5py"
    ]
)
