import mailbox
import os
from unicodedata import name
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
sig_event = pd.DataFrame(columns=['hi','Ra','Rm','et','mb','mj','Nb','HT','Rb','mh','rhb','rht','rlep','rhbh', 'ap', 'dp', 'mp'])
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
    #sig_event['pr'] = file["events"]["spanet_probability"][:]
    sig_event['ap'] = file["events"]["assignment_probability"][:]
    sig_event['dp'] = file["events"]["detection_probability"][:]
    sig_event['mp'] = file["events"]["marginal_probability"][:]


    flags = file["events"]["flags"][:]

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

num_trees = 500

#spliting trainning & evaluating
trainningcut = round(len(sig_event)*0.7)
sig_event_train = sig_event[:trainningcut]
sig_event_evaluate = sig_event[(trainningcut+1):]
flags_train = flags[:trainningcut]
flags_evaluate = flags[(trainningcut+1):]
# trainning 
dmatrix0 = xgb.DMatrix(data=sig_event_train, label=flags_train) # trainning
dmatrix1  = xgb.DMatrix(data=sig_event_evaluate, label=flags_evaluate) #testing

booster = xgb.train(xparams, dmatrix0, num_boost_round=num_trees)
trainingpredictionsx = booster.predict(dmatrix1)
print(trainingpredictionsx[1:10])

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Classifier Output")
for i in range(len(flags_evaluate)):
    if flags_train[i] == 0:
        bkg_pred.append(trainingpredictionsx[i])
    elif flags_train[i] == 1:
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
plt.savefig("classifieroutput_abs_spanet.png")

plt.figure()
plt.title("Classifier Output: >=6jets & >=4bjets")    
ntot=sum(n)
mtot=sum(m)

metric=wasserstein_distance(sig_pred, bkg_pred)
nai=833.9*0.288*0.001343*1000*300
mai=0.5085*0.288*0.1894*1000*300
#plt.yscale("log")
plt.hist(bincentres, bins=bins, weights=n*nai/ntot,alpha=0.5, label='ttbar SPANet, L=300 fb-1, integral='+str(sum(n)*nai/ntot))
plt.hist(bincentres, bins=bins, weights=m*mai/mtot,alpha=0.5, label='ttH SPANet, metric = '+str(round(metric,4))+'integral='+str(sum(m)*mai/mtot))

np.save('sig_bin',m*mai/mtot)
np.save('bkg_bin',n*nai/ntot)

plt.errorbar(bincentres, n/ntot, np.sqrt(n)/ntot, fmt='none', ecolor='b')
plt.errorbar(bincentres, m/mtot, np.sqrt(m)/mtot, fmt='none', ecolor='orange')

plt.legend()
plt.savefig("classifieroutput_normalized_spanet_tr_te.png")

plt.figure()
plt.cla()
auc = roc_auc_score(flags_evaluate, trainingpredictionsx)

print(round(auc,3))
    
fprx, tprx, _ = roc_curve(flags_evaluate.ravel(), trainingpredictionsx)
    
plt.plot(tprx, 1-fprx, label='XGBoost, AUC = '+str(round(auc,3)))
plt.title("ROC curves from SPANet BDT")
plt.xlabel('Signal Efficiency')
plt.ylabel('Background Rejection')
plt.legend()

plt.savefig('roc_spanet.png')

np.save('flag_span_w_sig_tr_te',flags_evaluate)
np.save('prob_span_w_sig_tr_te',trainingpredictionsx)
plt.cla()

xgb.plot_importance(booster, grid=False, importance_type='weight', title="Feature Importance: Weight >=6 jets(n_sig)")
plt.savefig('importance_weight.png')
plt.cla()
xgb.plot_importance(booster, grid=False, importance_type='gain', title="Feature Importance: Gain >=6 jets(n_sig)")
plt.savefig('importance_gain.png')
plt.cla()
xgb.plot_importance(booster, grid=False, importance_type='cover', title="Feature Importance: Cover >=6 jets(n_sig)")
plt.savefig('importance_cover.png')
plt.cla()
xgb.plot_importance(booster, grid=False, importance_type='total_gain', title="Feature Importance: Total Gain >=6 jets(n_sig)")
plt.savefig('importance_total_gain.png')
plt.cla()
xgb.plot_importance(booster, grid=False, importance_type='total_cover', title="Feature Importance: Total Cover >=6 jets(n_sig)")
plt.savefig('importance_total_cover.png')
plt.cla()
plt.hist(bincentres, bins=bins, weights=n*nai/ntot,alpha=0.5, label='ttbar SPANet, L=300 fb-1')
plt.savefig('ttbar_normalized.png')
plt.cla()
plt.hist(bincentres, bins=bins, weights=m*mai/mtot,alpha=0.5, label='ttH SPANet, L=300 fb-1')
plt.savefig('ttH_normalized.png')
# Save the model
booster.save_model('trainspanet0.xgb')
