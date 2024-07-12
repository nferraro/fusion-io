
import os, sys, glob
import numpy as np
from netCDF4 import Dataset
import pickle
import m3dc1_to_vmec_sfincs as convert

sys.path.append(os.environ['FIO_ROOT']+"util/sfincs")

m3dc1_outfile = './C1.h5'
te_start = 16000 # eV
te_end = 9000 # eV
timeslice = -1
ntheta = 400
nphi = 80
nsurf = 5 # Value is doubled in write_neo_input to account for the half mesh needed by VMEC
nfp = 2
ntor = 6
mpol = 8

convert.m3dc1_to_vmec(m3dc1_outfile, timeslice, nfp, ntheta, nphi, nsurf, vmec_mpol=mpol, vmec_ntor=ntor, te_start=te_start, te_end=te_end, dR0=3.24, vmec_file_descriptor="test", write_sfincs_input=True)
