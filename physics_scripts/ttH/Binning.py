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

