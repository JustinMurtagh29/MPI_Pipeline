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
import pickle
import pystan
import numpy as np
from os import path

if len(sys.argv) < 2:
    print("Data file input argument missing", file=sys.stderr)
    sys.exit(1)

in_file = sys.argv[1]
out_file = path.splitext(in_file)[0] + '.out.h5'

with h5py.File(in_file, 'r') as f:
    log10Asi1 = f['log10Asi1'][:].tolist()[0]
    log10Asi2 = f['log10Asi2'][:].tolist()[0]

if len(log10Asi1) != len(log10Asi2):
    print("ASI area vectors don't match up", file=sys.stderr)
    sys.exit(1)

if len(log10Asi1) == 0:
    print("ASI area vectors are empty", file=sys.stderr)
    sys.exit(1)

this_dir = path.dirname(path.realpath(__file__))
stan_file = path.join(this_dir, 'clusterConnections.stan')
pickle_file = path.join(this_dir, 'clusterConnections.pkl')

try:
    sm = pickle.load(open(pickle_file, 'rb'))
except:
    sm = pystan.StanModel(file=stan_file)
    pickle.dump(sm, open(pickle_file, 'wb'))

data = {'K': 3, 'N': len(log10Asi1), 'log10Asi1': log10Asi1, 'log10Asi2': log10Asi2}
fit = sm.sampling(data=data, seed=0, iter=1000, chains=1)
samples = fit.extract()

with h5py.File(out_file, 'w') as f:
    for key, data in samples.items():
        f.create_dataset(key, data=data)
