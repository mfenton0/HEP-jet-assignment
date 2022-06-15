import os
import numpy as np
import h5py
import xgboost as xgb
import pandas as pd
import matplotlib.pyplot as plt
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

sig_event = pd.DataFrame(columns=['hi','Ra','Rm','et','mb','mj','Nb','HT','Rb','mh','rhb','rht','rlep','rhbh'])
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

    flags = file["events"]["flags"][:]

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("higgs mass")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['hi'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['hi'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100
xrange = [0.0, 200]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("higgs mass")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("higgsmass.png")


sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['Ra'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['Ra'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 10]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("Average deltaR for all bjet pairs")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("Ravgbb.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['Rm'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['Rm'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 10]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("deltaR between bjet pairs with largest vector sumPt")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("Rmaxptbb.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['et'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['et'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 6]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("Max deltaEta between any 2 jets")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("etamaxjj.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['mb'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['mb'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 400]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("Mass of b jet pairs with smallest dR")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("massminRbb.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['mj'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['mj'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 400]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("Mass of jet pairs with smallest dR")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("massminRjj.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['Nb'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['Nb'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=12

xrange = [0.0, 12]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("Num of bjet pairs within 30GeV of H mass")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("Nbbhiggs30.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['HT'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['HT'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 3000]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("Scalar sum of Pt")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("HThad.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['Rb'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['Rb'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 10]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("dR between lepton and 2 bjets with smallest dR")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("Rminlbb.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['mh'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['mh'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 1000]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("mass of H and blep")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("mhblep.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['rhb'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['rhb'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 10]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("dR between 2bjet of H candidate")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("Rhbb.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['rht'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['rht'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 10]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("dR between H & ttbar")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("Rhtt.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['rlep'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['rlep'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 10]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("dR between H & tlep")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("Rhtlep.png")

sig_pred, bkg_pred = [], []
plt.figure()
plt.title("Average deltaR for all bjet pairs")
for i in range(len(flags)):
    if flags[i] == 0:
        bkg_pred.append(sig_event['rhbh'][i])
    elif flags[i] == 1:
        sig_pred.append(sig_event['rhbh'][i])
    else:
        print("ERROR, CHECK YOUR FLAGS")

nbins=100

xrange = [0.0, 10]
binwidth = (xrange[1]-xrange[0])/nbins
norm = False
errors=True
# the first plot we specify the binning and get them from the return for the rest
n, bins,patches = plt.hist(bkg_pred,alpha=0.5,label='ttbar',bins=nbins,range=xrange, density=norm)
bincentres = bins[:-1]+0.5*binwidth # bins[-1] since bin length is nbins+1
m, bins, _ = plt.hist(sig_pred, alpha=0.5, label="ttH", bins=bins, density=norm)

plt.figure()
plt.title("dR between H & bhad")    
ntot=sum(n)
mtot=sum(m)

plt.hist(bincentres, bins=bins, weights=n/ntot,alpha=0.5, label='ttbar')
plt.hist(bincentres, bins=bins, weights=m/mtot,alpha=0.5, label='ttH')

plt.legend()
plt.savefig("Rhbhad.png")