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
sig_event = pd.DataFrame(columns=['hi','Ra','Rm','et','mb','mj','Nb','HT','Rb','mh','rhb','rht','rlep','rhbh', 'ap', 'dp', 'mp', 'pr'])
with h5py.File("spanet_features.h5",'r') as file:
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
    sig_event['pr'] = file["events"]["spanet_probability"][:]
    sig_event['ap'] = file["events"]["assignment_probability"][:]
    sig_event['dp'] = file["events"]["detection_probability"][:]
    sig_event['mp'] = file["events"]["marginal_probability"][:]


    flags = file["events"]["flags"][:]

# Set the xgboost hyperparameters
# https://xgboost.readthedocs.io/en/latest/parameter.html


sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Classifier Output")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['pr'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['pr'][i])
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
plt.savefig("classifieroutput_abs_spanet.png")

plt.figure()
plt.title("Classifier Output: >=6jets & >=4bjets")    
ntot=sum(n)
mtot=sum(m)

metric=wasserstein_distance(sig_pred, bkg_pred)
plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar SPANet')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH SPANet, metric = '+str(round(metric,4)))

plt.errorbar(bincentres, n/ntot, np.sqrt(n)/ntot, fmt='none', ecolor='b')
plt.errorbar(bincentres, m/mtot, np.sqrt(m)/mtot, fmt='none', ecolor='orange')

plt.legend()
plt.savefig("classifieroutput_normalized_spanet.png")

plt.figure()
plt.cla()
auc = roc_auc_score(flags, sig_event['pr'])
    
fprx, tprx, _ = roc_curve(flags.ravel(), sig_event['pr'])
    
plt.plot(tprx, 1-fprx, label='SPANet only, AUC = '+str(round(auc,3)))
plt.title("ROC curves from SPANet >=6jets & >=4bjets")
plt.xlabel('Signal Efficiency')
plt.ylabel('Background Rejection')
plt.legend()

flags1 = np.load('flag_span_w_sig.npy')
prob1 = np.load('prob_span_w_sig.npy')
print(flags1)
auc1 = roc_auc_score(flags1, prob1)
    
fprx, tprx, _ = roc_curve(flags1.ravel(), prob1)
    
plt.plot(tprx, 1-fprx, label='SPANet BDT(w_sig), AUC = '+str(round(auc1,3)))
plt.legend()

flags2 = np.load('flag_span_n_sig.npy')
prob2 = np.load('prob_span_n_sig.npy')
print(flags1)
auc2 = roc_auc_score(flags2, prob2)
    
fprx, tprx, _ = roc_curve(flags2.ravel(), prob2)
    
plt.plot(tprx, 1-fprx, label='SPANet BDT(n_sig), AUC = '+str(round(auc2,3)))
plt.legend()
plt.savefig('roc_spanet.png')

flags3 = np.load('flag_klf_sig.npy')
prob3 = np.load('prob_klf_sig.npy')

auc3 = roc_auc_score(flags3, prob3)
    
fprx, tprx, _ = roc_curve(flags3.ravel(), prob3)
    
plt.plot(tprx, 1-fprx, label='KLFitter, AUC = '+str(round(auc3,3)))
plt.legend()
plt.savefig('roc_spanet.png')

# Save the model

