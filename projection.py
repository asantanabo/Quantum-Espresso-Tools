#!/usr/bin/env python

import numpy as np
import sys,os,re
import itertools
import matplotlib.pyplot as plt

# Functions to be used by the program

def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected.
        Required arguments:                                                                                                                                                                                     
        l -- The iterable to be sorted.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

def chunker(iterable, chunksize):
    return list(map(None,*[iter(iterable)]*chunksize))

def read_qpoint(filename):
    freqs = []
    disps = []
    with open(filename) as f:

        for line in f:

            if line.strip():
                data = line.strip().split()
            else:
                continue

            if "THz" in line:
                freq = float(data[3])
                freqs.append(freq)
            
            else:
                data = list(map(float, data))
                disps.append(data)

    freqs = np.array(freqs)
    nfreqs = len(freqs)
    disps = np.array(disps).reshape(nfreqs, -1, 6)

    return freqs, disps

#### END-FUNCTIONS  ################################

# Processing the q-point files

# try:
#     infilename = sys.argv[1];
# except:
#     print "Usage",sys.argv[0], "outfile_from_quantum_espresso_eigenvectors"; sys.exit(1)
# 
# ifile = open(infilename, 'r')
# test = [];

# Organizing the q-points-i.dat files
current = os.getcwd()
files = []
for file in os.listdir(current):
    if file.endswith(".dat"):
        files.append(os.path.join(current, file))

files = sorted_nicely(files)

allfreqs = []
alldisps = []
for qfile in files:
    freqs, disps = read_qpoint(qfile)
    allfreqs.append(freqs)
    alldisps.append(disps)

allfreqs = np.array(allfreqs)
nfreqs = len(allfreqs[0])

disp0 = np.array(alldisps[0])
alldisps = np.array(alldisps[1:])

result = []
for mode in disp0:
    for qpoint in alldisps:
        for mode1 in qpoint:
            A = (mode * mode1).sum()
            result.append(A)

result = np.array(result).reshape(-1, nfreqs, nfreqs)
amp_matrix = result.reshape((144,144))

plt.matshow(amp_matrix, vmin=0, vmax=1.0, cmap='hot')
plt.colorbar()
plt.xlim(0,20)
plt.ylim(20,0)
plt.grid(True)
plt.show()

