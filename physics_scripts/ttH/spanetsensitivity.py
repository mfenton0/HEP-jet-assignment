from pathlib import Path
from math import sqrt

import numpy as np
import matplotlib.pyplot as plt

import pyhf
import pyhf.readxml
from pyhf.contrib.viz import brazil

bkg11 = np.load('span_bkg_bin.npy')
sig11 = np.load('span_sig_bin.npy')
ntot=sum(bkg11)
mtot = sum(sig11)
nai=833.9*0.288*0.001343*1000*300
mai=0.5085*0.288*0.1894*1000*300
uncer = []
sig1 = []
bkg1 = []
data1 = []
for i in range(0,len(bkg11)):
    uncer.append(0.2*sqrt(bkg11[i])*nai/ntot)
    sig1.append(sig11[i]*mai/mtot)
    bkg1.append(bkg11[i]*nai/ntot) 
    data1.append(bkg11[i]*nai/ntot)

pyhf.set_backend("numpy")
model = pyhf.simplemodels.uncorrelated_background(
    signal=sig1, bkg=bkg1, bkg_uncertainty= uncer
)
data = data1 + model.config.auxdata

poi_vals = np.linspace(0, 0.6, 41)


(
    obs_limit,
    exp_limits,
    (poi_tests, tests),
) = pyhf.infer.intervals.upper_limits.upper_limit(
    data, model, poi_vals, level=0.05, return_results=True
)

fig, ax = plt.subplots()
fig.set_size_inches(7, 5)
brazil.plot_results(poi_vals, tests, test_size=0.05, ax=ax)
fig.show()
plt.title("SPANet BDT")
fig.savefig("sense.png")
