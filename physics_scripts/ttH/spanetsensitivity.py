from pathlib import Path
from math import sqrt

import numpy as np
import matplotlib.pyplot as plt

import pyhf
import pyhf.readxml
from pyhf.contrib.viz import brazil

bkg11 = np.load('bkg_bin.npy')
sig11 = np.load('sig_bin.npy')
uncer = []
sig1 = []
bkg1 = []
data1 = []
for i in range(0,len(bkg11)):
    uncer.append(sqrt(bkg11[i]))
    sig1.append(sig11[i])
    bkg1.append(bkg11[i]) 
    data1.append(bkg11[i])

pyhf.set_backend("numpy")
model = pyhf.simplemodels.uncorrelated_background(
    signal=sig1, bkg=bkg1, bkg_uncertainty= uncer
)
data = data1 + model.config.auxdata

poi_vals = np.linspace(0, 0.6, 41)
results = [
    pyhf.infer.hypotest(
        test_poi, data, model, test_stat="qtilde", return_expected_set=True
    )
    for test_poi in poi_vals
]

fig, ax = plt.subplots()
fig.set_size_inches(7, 5)
brazil.plot_results(poi_vals, results, ax=ax)
fig.show()
plt.title("SPANet BDT")
fig.savefig("sense.png")
