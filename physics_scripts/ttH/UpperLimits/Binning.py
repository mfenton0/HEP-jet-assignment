import os
import numpy as np
import h5py
import xgboost as xgb
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wasserstein_distance
from sklearn.metrics import roc_curve, roc_auc_score

def Autobinning(flagfile,probfile,nbins=5):
    flag = np.load(flagfile)
    prob = np.load(probfile)

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
    #Za = nbins*(Nsig/len(flag))
    #Zb = nbins-Za
    Za = 0.5*nbins
    Zb = 0.5*nbins
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
            #print(BDTcut)

    Interval.append(BDTcut)
    Interval.append(0)
    print(flagfile)
    print(Interval)

Autobinning('flag_klf_sig_tr_te.npy','prob_klf_sig_tr_te.npy',5)
Autobinning('flag_dnn_w_sig_tr_te.npy','prob_dnn_w_sig_tr_te.npy',5)
Autobinning('flag_span_w_sig_tr_te_BDT.npy','prob_span_w_sig_tr_te_BDT.npy',5)
Autobinning('flag_spanet_pretraining.npy','prob_spanet_pretraining.npy',5)
Autobinning('flag_spanet_finetuning.npy','prob_spanet_finetuning.npy',5)

Autobinning('flag_klf_sig_tr_te.npy','prob_klf_sig_tr_te.npy',10)
Autobinning('flag_dnn_w_sig_tr_te.npy','prob_dnn_w_sig_tr_te.npy',10)
Autobinning('flag_span_w_sig_tr_te_BDT.npy','prob_span_w_sig_tr_te_BDT.npy',10)
Autobinning('flag_spanet_pretraining.npy','prob_spanet_pretraining.npy',10)
Autobinning('flag_spanet_finetuning.npy','prob_spanet_finetuning.npy',10)

Autobinning('flag_klf_sig_tr_te.npy','prob_klf_sig_tr_te.npy',8)
Autobinning('flag_dnn_w_sig_tr_te.npy','prob_dnn_w_sig_tr_te.npy',8)
Autobinning('flag_span_w_sig_tr_te_BDT.npy','prob_span_w_sig_tr_te_BDT.npy',8)
Autobinning('flag_spanet_pretraining.npy','prob_spanet_pretraining.npy',8)
Autobinning('flag_spanet_finetuning.npy','prob_spanet_finetuning.npy',8)