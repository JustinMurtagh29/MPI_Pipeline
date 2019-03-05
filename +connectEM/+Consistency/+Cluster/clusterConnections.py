#!/usr/bin/env python3

# This Python script serves as shim between MATLAB and Stan. It takes
# the path to a HDF5 file, loads the log10(axon-spine interface areas)
# and passes them to Stan to fit the connection type parameters.
#
# For more information about the connection type model, see the
# clusterConnections.stan file in the same directory.
#
# Written by
#   Alessandro Motta <alessandro.motta@brain.mpg.de>

import sys
import h5py
import pystan
import numpy as np
from os import path

if len(sys.argv) < 2:
    print("Data file input argument missing", file=sys.stderr)
    sys.exit(1)

f = h5py.File(sys.argv[1], 'r')
log10Asi1 = f['log10Asi1'][:].tolist()[0]
log10Asi2 = f['log10Asi2'][:].tolist()[0]
f.close()

if len(log10Asi1) != len(log10Asi2):
    print("ASI area vectors don't match up", file=sys.stderr)
    sys.exit(1)

if len(log10Asi1) == 0:
    print("ASI area vectors are empty", file=sys.stderr)
    sys.exit(1)

this_dir = path.dirname(path.realpath(__file__))
sm = pystan.StanModel(file=path.join(this_dir, "clusterConnections.stan"))
# TODO: Pickle compiled model

# TODO: Bump number of iterations
data = {'K': 3, 'N': len(log10Asi1), 'log10Asi1': log10Asi1, 'log10Asi2': log10Asi2}
fit = sm.sampling(data=data, iter=1000, chains=4)
print(fit)
