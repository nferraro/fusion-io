1
import os, sys, glob
import numpy as np
from netCDF4 import Dataset
import pickle
#import m3dc1_to_vmec_sfincs as convert

sys.path.append(os.environ['FIO_ROOT']+"/util/vmec-sfincs")
import m3dc1_to_vmec_sfincs as convert


#~/src/fusion-io/util/_stellar/write_neo_input -bootstrap 1 -m3dc1 ../C1.h5 0 -te_start 9400 -te_end 850 -nr 40 -nphi 720 -ntheta 400 -tol 1 -R0 12.5
m3dc1_outfile = '../../C1.h5'
te_start =  9400
te_end = 1200
usingT = 1
bootstrap = 1
timeslice = 0
ntheta = 400
nphi = 540
nsurf =20 # Value is doubled in write_neo_input to account for the half mesh needed by VMEC
nfp = 2
ntor = 12
mpol = 16
R =12
tol=1
psi_t_end_input=1#if 0 use psi_t_end of neo_input.nc that was created

if(usingT == 1):
    convert.m3dc1_to_vmec(m3dc1_outfile, timeslice, nfp, ntheta, nphi, nsurf, vmec_mpol=mpol, vmec_ntor=ntor, bootstrap=bootstrap,te_start=te_start, te_end=te_end, R0=R, tol=tol, usingT=usingT,vmec_file_descriptor="test", write_sfincs_input=True,psi_t_end_input=psi_t_end_input)
else:
    convert.m3dc1_to_vmec(m3dc1_outfile, timeslice, nfp, ntheta, nphi, nsurf, vmec_mpol=mpol, vmec_ntor=ntor, bootstrap=bootstrap, psi_start=psi_start, psi_end=psi_end, tol=tol, dR0=dR, usingT=usingT,vmec_file_descriptor="test", write_sfincs_input=True,psi_t_end_input=psi_t_end_input)
