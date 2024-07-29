
import os, sys, glob
import numpy as np
from netCDF4 import Dataset
import pickle
#import m3dc1_to_vmec_sfincs as convert

sys.path.append(os.environ['FIO_ROOT']+"/util/vmec-sfincs")
import m3dc1_to_vmec_sfincs as convert

#~/src/fusion-io/util/_stellar/write_neo_input -bootstrap 1 -m3dc1 ../C1.h5 1 -te_start 2413 -te_end 1000 -nr 25 -nphi 100 -dR0 1.8
#-bootstrap 1 -m3dc1 ../C1.h5 1 -te_start 2413 -te_end 1000 -nr 25 -nphi 100 -dR0 1.8
m3dc1_outfile = '../C1.h5'
#te_start =  2413
#te_end = 1000
usingT = 0
bootstrap = 0
psi_start =  0.7
psi_end = 0.96 
timeslice = 0
ntheta = 400
nphi = 100
nsurf = 100 # Value is doubled in write_neo_input to account for the half mesh needed by VMEC
nfp = 2
ntor = 6
mpol = 8
dR =0
tol=0.001

if(usingT == 1):
    convert.m3dc1_to_vmec(m3dc1_outfile, timeslice, nfp, ntheta, nphi, nsurf, vmec_mpol=mpol, vmec_ntor=ntor, bootstrap=bootstrap,te_start=te_start, te_end=te_end, dR0=dR, tol=tol, vmec_file_descriptor="test", write_sfincs_input=True)
else:
    convert.m3dc1_to_vmec(m3dc1_outfile, timeslice, nfp, ntheta, nphi, nsurf, vmec_mpol=mpol, vmec_ntor=ntor, bootstrap=bootstrap, psi_start=psi_start, psi_end=psi_end, tol=tol, dR0=dR, vmec_file_descriptor="test", write_sfincs_input=True)
