import os
import numpy as np
import h5py
import xgboost as xgb
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wasserstein_distance
from sklearn.metrics import roc_curve, roc_auc_score

DATA_STRUCTURE = np.dtype(
    [
        ("higgsm", "f4"),
        ("Ravgbb", "f4"),
        ("Rmaxptbb", "f4"),
        ("etamaxjj", "f4"),
        ("massminRbb", "f4"),
        ("massminRjj", "f4"),
        ("Nbbhiggs30", "f4"),
        ("HThad", "f4"),
        ("Rminlbb", "f4"),
        ("mhblep", "f4"),
        ("Rhbb", "f4"),
        ("Rhtt", "f4"),
        ("Rhtlep", "f4"),
        ("Rhbhad", "f4"),
    ]
)
sig_event = pd.DataFrame(columns=['hi','Ra','Rm','et','mb','mj','Nb','HT','Rb','mh','rhb','rht','rlep','rhbh', 'li'])

with h5py.File("output_signal.h5",'r') as file:
    data_sig = file["events"]["Nbbhiggs30"][:]
    
    sig_event['hi'] = file["events"]["higgsm"][:]
    sig_event['Ra'] = file["events"]["Ravgbb"][:]
    sig_event['Rm'] = file["events"]["Rmaxptbb"][:]
    sig_event['et'] = file["events"]["etamaxjj"][:]
    sig_event['mb'] = file["events"]["massminRbb"][:]
    sig_event['mj'] = file["events"]["massminRjj"][:]
    sig_event['Nb'] = file["events"]["Nbbhiggs30"][:]
    sig_event['HT'] = file["events"]["HThad"][:]
    sig_event['Rb'] = file["events"]["Rminlbb"][:]
    sig_event['mh'] = file["events"]["mhblep"][:]
    sig_event['rhb'] = file["events"]["Rhbb"][:]
    sig_event['rht'] = file["events"]["Rhtt"][:]
    sig_event['rlep'] = file["events"]["Rhtlep"][:]
    sig_event['rhbh'] = file["events"]["Rhbhad"][:]
    sig_event['li'] = file["events"]["likelihoodmax"][:]

    flags = file["events"]["spanet_signal_target"][:]



flag = []
for j in range(len(flags)):
    if flags[j] == 0:
        flag.append(0)
    else:
        flag.append(1)

b = []
for jj in range(len(flag)):
    b.append(np.array(flag[jj]))
flag = np.array(b)
# Set the xgboost hyperparameters
# https://xgboost.readthedocs.io/en/latest/parameter.html
xparams = {}

xparams['learning_rate'] = 0.01 
xparams['max_depth'] = 3
xparams['objective'] = 'binary:logistic'
xparams['seed'] = 42
xparams['eval_metric'] = 'auc'
#xparams['nthread'] = 8
xparams=list(xparams.items())

num_trees = 7000
trainningcut = round(len(sig_event)*0.7)
sig_event_train = sig_event[:trainningcut]
sig_event_evaluate = sig_event[(trainningcut+1):]
flag_train = flag[:trainningcut]
flag_evaluate = flag[(trainningcut+1):]
# trainning 
dmatrix0 = xgb.DMatrix(data=sig_event_train, label=flag_train) # trainning
dmatrix1  = xgb.DMatrix(data=sig_event_evaluate, label=flag_evaluate) #testing

booster = xgb.train(xparams, dmatrix0, num_boost_round=num_trees)
trainingpredictionsx = booster.predict(dmatrix1)
print(trainingpredictionsx[1:10])

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Classifier Output")
for i in range(len(flag_evaluate)):
    if flag_evaluate[i] == 0:
        bkg_pred.append(trainingpredictionsx[i])
    elif flag_evaluate[i] == 1:
        sig_pred.append(trainingpredictionsx[i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=20
xrange = [0.0, 1.0]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar Test',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
if errors:
    plt.errorbar(bincentres,n,np.sqrt(n), fmt='none', ecolor='b')
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH Test", bins=bins, density=norm)

if errors:
    plt.errorbar(bincentres, m,np.sqrt(m), fmt='none', ecolor='orange')

plt.legend()
plt.savefig("classifieroutput_abs.png")

plt.figure()
plt.title("Classifier Output: ==8jets & >=4bjets")    
ntot=sum(n)
mtot=sum(m)

metric=wasserstein_distance(sig_pred, bkg_pred)
plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar KLFitter')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH KLFitter, metric = '+str(round(metric,4)))

plt.errorbar(bincentres, n/ntot, np.sqrt(n)/ntot, fmt='none', ecolor='b')
plt.errorbar(bincentres, m/mtot, np.sqrt(m)/mtot, fmt='none', ecolor='orange')

plt.legend()
plt.savefig("classifieroutput_normalized.png")

plt.figure()
plt.cla()
plt.title("ROC curves from KLFitter BDT")   
auc = roc_auc_score(flag_evaluate, trainingpredictionsx)
    
fprx, tprx, _ = roc_curve(flag_evaluate.ravel(), trainingpredictionsx)
    
plt.plot(tprx, 1-fprx, label='XGBoost, AUC = '+str(round(auc,3)))
plt.xlabel('Signal Efficiency')
plt.ylabel('Background Rejection')
plt.legend()

plt.savefig('roc.png')

np.save('flag_klf_sig',flag)
np.save('prob_klf_sig',trainingpredictionsx)
# Save the model
booster.save_model('trainon0.xgb')
