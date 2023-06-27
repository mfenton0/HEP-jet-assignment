from pathlib import Path
from math import sqrt

import numpy as np
import matplotlib.pyplot as plt

import pyhf
import pyhf.readxml
from pyhf.contrib.viz import brazil

SMxec = 0.507
unit = 'pb'
nbins=5
bkg11 = np.load('dnn_bkg_bin_nbins'+str(nbins)+'.npy')
sig11 = np.load('dnn_sig_bin_nbins'+str(nbins)+'.npy')
BkgError = sqrt(490*490+210*210+130*130+61*61+27*27+49*49+24*24)
kttbb = 1
ntot=sum(bkg11)
mtot = sum(sig11)
lum = 139
#nai=833.9*0.288*0.001343*1000*lum
#mai=0.507*0.288*0.577*0.1894*1000*lum
nai = 3160 + 1530 + 720 + 215 + 55 + 246 + 55
mai = 213 + 113 + 59.9 + 13.9 + 3.09 + 35.1 +8.5
BkgError = sqrt(810*810+360*360+180*180+72*72+29*29+120*120+38*38)
nai = 4350+2100+1000+301+80+470+117
bkgrate = BkgError/nai
mai = 502
nai = 4337
kttbb = 1.5
if nbins==5:
    bkgrate = 0.5
    mai = mai*140/139
    nai = nai*140/139
else:
    bkgrate = 0.3
    mai = mai*300/139
    nai = nai*300/139
uncer = []
sig1 = []
bkg1 = []
data1 = []
alpha = 0.05
for i in range(0,len(bkg11)):
    uncer.append(bkgrate*sqrt(bkg11[i])*nai/ntot)
    sig1.append(sig11[i]*mai/mtot)
    bkg1.append(kttbb*bkg11[i]*nai/ntot) 
    data1.append(kttbb*bkg11[i]*nai/ntot+sig11[i]*mai/mtot)
pyhf.set_backend("numpy")
model = pyhf.simplemodels.uncorrelated_background(
    signal=sig1, bkg=bkg1, bkg_uncertainty= uncer
)
data = data1 + model.config.auxdata

poi_vals = np.linspace(0, 2, 80)


(
    obs_limit,
    exp_limits,
    (poi_tests, tests),
) = pyhf.infer.intervals.upper_limits.upper_limit(
    data, model, poi_vals, level=alpha, return_results=True
)

fig, ax = plt.subplots()
fig.set_size_inches(7, 5)
brazil.plot_results(poi_vals, tests, test_size=alpha, ax=ax)
fig.show()
plt.title("PermutationDNN BDT, expected limits:"+str(np.round(exp_limits[2],3))+" bkg uncertainty: 20% (300fb-1)")
fig.savefig("sensednn.png")

asymp_calc = pyhf.infer.calculators.AsymptoticCalculator(
    data, model, test_stat='qtilde'
)
teststat = asymp_calc.teststatistic(poi_test=obs_limit)
print("DNNBDT:"+str(np.round(exp_limits[2],3))+" xsec:"+str(np.round((np.round(exp_limits[2],3)*SMxec),3))+"(pb) "+ f'qtilde = {teststat}')
