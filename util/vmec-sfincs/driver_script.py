
import os, sys, glob
import numpy as np
from netCDF4 import Dataset
import pickle
import m3dc1_to_vmec_sfincs as convert

sys.path.append(os.environ['FIO_ROOT']+"util/sfincs")

m3dc1_outfile = '../C1.h5'
te_start = 1200 # eV
te_end = 100 # eV
timeslice = -1
ntheta = 400
nphi = 32
nsurf = 50 # Value is doubled in write_neo_input to account for the half mesh needed by VMEC
nfp = 1
ntor = 6
mpol = 8
dR0 = 0.83

convert.m3dc1_to_vmec(m3dc1_outfile, timeslice, nfp, ntheta, nphi, nsurf, vmec_mpol=mpol, vmec_ntor=ntor, te_start=te_start, te_end=te_end, dR0=dR0, vmec_file_descriptor="new_tcv", write_sfincs_input=True)
