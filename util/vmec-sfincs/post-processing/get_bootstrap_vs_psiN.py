import os, sys, glob
import math
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import h5py

# This script takes the output of up to 100 SFINCS runs, and outputs the bootstrap current from each, in SI units

# Normalizations
Bbar = 1.0 # T
nbar = 1.e20 # m^-3
vbar = 4.376947085e+05 # m/s
e = 1.602176565e-19


# The format for all directories should be of the form "XX" where "XX" is a numeral between 00-99
sfincs_runs = glob.glob("./[0-9][0-9]")
nruns = len(sfincs_runs)

psiN = []
FSABJ_SI = []

outputFilename = "sfincsOutput.h5"
for irun in range(nruns):
    filename = sfincs_runs[irun]+"/"+outputFilename
    print("irun = {}".format(sfincs_runs[irun]))
    if not os.path.isfile(filename):
        print("Directory "+sfincs_runs[irun]+" does not have a "+outputFilename+" file. Exiting...")
        exit()

    try:
        f = h5py.File(filename, 'r')
    except IOError:
        print("Unable to open "+filename+" even though this file exists. Skipping...")
        continue

    try:
        finished = f["finished"]
    except KeyError:
        print("Run in directory"+sfincs_runs[irun]+" does not appear to have finished. Skipping...")
        continue
        
    psiN.append(f["psiN"][()])
    FSABJ_SI.append(f["FSABjHat"][()][0] * Bbar * nbar * vbar * e * 1.e-6) # Convert to MA

results = sorted(zip(psiN, FSABJ_SI), key = lambda x: x[0])
psiN_final = [x for x,_ in results]
FSABJ_final = [x for _,x in results]


outfile = open("./FSABJ.dat", 'w')
outfile.write("Row 1: psiN, row 2: <B.J> (MA T/m^2)\n")
for iline in range(len(psiN_final)):
    outfile.write(str(psiN_final[iline])+" ")
outfile.write("\n")
for iline in range(len(psiN_final)):
    outfile.write(str(FSABJ_final[iline])+" ")
outfile.write("\n")
outfile.close()

plt.figure()
plt.plot(np.sqrt(psiN_final), FSABJ_final)
plt.show()
exit()
