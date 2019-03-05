import h5py
import numpy as np
import pystan

f = open('/home/amotta/Desktop/gaussian_mixture.stan', 'r')
model_code = f.read()
f.close()

f = h5py.File('/home/amotta/Desktop/gaussian_mixture_data.h5', 'r')
log10Asi1 = f['log10Area1'][:].tolist()[0]
log10Asi2 = f['log10Area2'][:].tolist()[0]
f.close()

sm = pystan.StanModel(model_code=model_code)
data = {'K': 3, 'N': 5290, 'log10Asi1': log10Asi1, 'log10Asi2': log10Asi2}
fit = sm.sampling(data=data, iter=1000, chains=1)
