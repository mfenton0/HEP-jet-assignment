import os
import numpy as np
import h5py
import xgboost as xgb
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wasserstein_distance
from sklearn.metrics import roc_curve, roc_auc_score

flag = np.load('flag_spanet_onlyfr.npy')
prob = np.load('prob_spanet_onlyfr.npy')

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Classifier Output")
for i in range(len(flag)):
    if flag[i] == 0:
        bkg_pred.append(prob[i])
    elif flag[i] == 1:
        sig_pred.append(prob[i])
    else:
        print("ERROR, CHECK YOUR FLAGS")
        
Nsig=len(sig_pred)
Nbkg=len(bkg_pred)
nbins=5
Za = nbins*(Nsig/len(flag))
Zb = nbins-Za
def Zfunction(nsig,nbkg,Nsig,Nbkg):
    return ((Za*nsig/Nsig)+(Zb*nbkg/Nbkg))



sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Classifier Output")
for i in range(len(flag)):
    if flag[i] == 0:
        bkg_pred.append(prob[i])
    elif flag[i] == 1:
        sig_pred.append(prob[i])
    else:
        print("ERROR, CHECK YOUR FLAGS")


xrange = [0.0, 1.0]
BDTcut = xrange[1]
BDTUppercut = xrange[1]

Interval = []
while len(Interval)<(nbins-1):
    Interval.append(BDTcut)
    BDTUppercut = Interval[-1]
    Zbin = 0
    while Zbin<1:
        BDTcut = BDTcut-0.01
        signum = 0
        bkgnum = 0
        for i in range(len(sig_pred)):
            if((sig_pred[i]>BDTcut)&(sig_pred[i]<BDTUppercut)):
                signum += 1
        for i in range(len(bkg_pred)):
            if((bkg_pred[i]>BDTcut)&(bkg_pred[i]<BDTUppercut)):
                bkgnum += 1
        Zbin = Zfunction(signum,bkgnum,Nsig,Nbkg)
        print(BDTcut)

    
    print("one bin formed")
Interval.append(BDTcut)
Interval.append(0)
print(Interval)
"""binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
ibins = np.array([0, 0.45, 0.5, 0.55, 0.6, 1])
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar KLFitter',bins=ibins, density=norm)
#bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
bincentres = np.array([0.225,0.457,0.525,0.575,0.8])
if errors:
    plt.errorbar(bincentres,n,np.sqrt(n), fmt='none', ecolor='b')
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH KLFitter",bins=ibins,density=norm)

if errors:
    plt.errorbar(bincentres, m,np.sqrt(m), fmt='none', ecolor='orange')

plt.legend()
plt.savefig("classifieroutput_abs.png")

plt.figure()
plt.title("Classifier Output: >=6jets & >=4bjets")    
ntot=sum(n)
mtot=sum(m)

metric=wasserstein_distance(sig_pred, bkg_pred)
nai=833.9*0.288*0.001343*1000*300
mai=0.5085*0.288*0.1894*1000*300
#plt.yscale("log")
#plt.hist(bincentres, bins=bins, weights=n*nai/ntot,alpha=0.5, label='ttbar KLFitter, L=300 fb-1, integral='+str(sum(n)*nai/ntot))
#plt.hist(bincentres, bins=bins, weights=m*mai/mtot,alpha=0.5, label='ttH KLFitter, metric = '+str(round(metric,4))+'integral='+str(sum(m)*mai/mtot))
plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar KLFitter BDT'+str(nbins)+'bins',histtype="step")
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH KLFitter BDT metric = '+str(round(metric,4)),histtype="step")
np.save('klf_sig_bin_nbins'+str(nbins),m)
np.save('klf_bkg_bin_nbins'+str(nbins),n)

plt.errorbar(bincentres, n/ntot, np.sqrt(n)/ntot, fmt='none', ecolor='b')
plt.errorbar(bincentres, m/mtot, np.sqrt(m)/mtot, fmt='none', ecolor='orange')

plt.legend()
plt.savefig("classifieroutput_normalized"+str(nbins)+".png")

# Save the model"""

