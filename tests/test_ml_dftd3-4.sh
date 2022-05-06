# Quick test for DFT-D3
# Note that the tests are just for checking the implementation and most likely provide nonsensical results.
# 1st step: install dftd3
# conda install dftd3 -c psi4
#Execute the following
# we replaced the calculator for ML+vdW with the d3 calculator, which is implemented in ase. It uses the ML model as an input for structure relaxation and the dftd3 correction on top
mkdir OptAuC
python test_ml_d3.py AuC_geometry.in OptAuC MLModels/Au_C/EF --qbc --nmodels 4 --fmax 0.1

# for dft 4 install dftd4 according to https://github.com/dftd4/dftd4
mkdir OptAu50
python test_ml.py Au50.in OptAu50 MLModels/Au_C/EF --fmax 0.1 --dftd4 --functional_d4 pbe
# note that this test is only done for Au50 cluster without diamond surface due to computational efficiency
