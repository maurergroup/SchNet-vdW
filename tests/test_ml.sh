# Quick test

#Au@C
#python test_ml.py AuC_geometry.in OptAuC MLModels/Au_C/EF/ --hirshfeld_modelpath MLModels/Au_C/Hirshfeld/ --vdw mbd --ts TS --nmodels 4 --qbc > OptAuC/opt.log

# To run bh, simply add "--mode bh". The default is "opt".
# e.g.:
# python test_ml.py AuC_geometry.in OptAuC MLModels/Au_C/EF --hirshfeld_modelpath MLModels/Au_C/Hirshfeld --vdw mbd --ts TS --nmodels 4 --qbc --mode bh > OptAuC/opt.log

#X2O

#Get Hessian file
#The Hessian input was created by Jürgen Wieferink and is available via the FHI-aims package (Lindh.py in utilities, python2). 
#A few other versions of the lindh script (including python3 version) is available via the package: https://github.com/sabia-group/gensec DOI 10.5281/zenodo.5651513. 

#run optimization
python test_ml.py B2O.in OptX2O MLModels/X2O/EF --hirshfeld_modelpath MLModels/X2O/Hirshfeld/ --vdw vdw --ts TSSurf --Hessianfile OptX2O/HessianB2O.dat --kgrid 7 7 1 > OptX2O/opt.log
#use more than 1 model:
#python test_ml.py B2O.in OptX2O MLModels/X2O/EF --hirshfeld_modelpath MLModels/X2O/Hirshfeld/ --vdw vdw --ts TSSurf --Hessianfile OptX2O/HessianB2O.dat --kgrid 7 7 1 --qbc --nmodels 4 > OptX2O/opt.log


