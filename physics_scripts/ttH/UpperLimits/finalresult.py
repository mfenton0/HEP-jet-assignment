import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import wasserstein_distance
from sklearn.metrics import roc_curve, roc_auc_score
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "monospace",
    "font.monospace": 'Computer Modern Typewriter'
})
import seaborn as sb

sb.set_theme("poster", "whitegrid", palette='deep', font_scale=1.1, rc={'figure.figsize':(12,10), 'text.usetex': True})
def getbincentre(bins):
    bincentre = []
    for i in range(0,len(bins)-1):
        bincentre.append(0.5*bins[i]+0.5*bins[i+1])
    return bincentre


flags1 = np.load('flag_klf_sig_tr_te.npy')
prob1 = np.load('prob_klf_sig_tr_te.npy')

klfsig = []
klfbkg = []
for i in range(len(flags1)):
    if flags1[i]==0:
        klfbkg.append(prob1[i])
    else:
        klfsig.append(prob1[i])

plt.figure()
plt.cla()

auc1 = roc_auc_score(flags1, prob1)
    
fprx, tprx, _ = roc_curve(flags1.ravel(), prob1)
    
plt.title(" ")
plt.xlabel('Signal Efficiency')
plt.ylabel('Background Rejection')

flags2 = np.load('flag_dnn_w_sig_tr_te.npy')
prob2 = np.load('prob_dnn_w_sig_tr_te.npy')

dnnsig = []
dnnbkg = []
for i in range(len(flags2)):
    if flags2[i]==0:
        dnnbkg.append(prob2[i])
    else:
        dnnsig.append(prob2[i])


auc2 = roc_auc_score(flags2, prob2)
    
fprx2, tprx2, _ = roc_curve(flags2.ravel(), prob2)

flags3 = np.load('flag_span_w_sig_tr_te_BDT.npy')
prob3 = np.load('prob_span_w_sig_tr_te_BDT.npy')

spansig = []
spanbkg = []
for i in range(len(flags3)):
    if flags3[i]==0:
        spanbkg.append(prob3[i])
    else:
        spansig.append(prob3[i])

auc3 = roc_auc_score(flags3, prob3)
    
fprx3, tprx3, _ = roc_curve(flags3.ravel(), prob3)

flags4 = np.load('flag_spanet_finetuning.npy')
prob4 = np.load('prob_spanet_finetuning.npy')

ftsig = []
ftbkg = []
for i in range(len(flags4)):
    if flags4[i]==0:
        ftbkg.append(prob4[i])
    else:
        ftsig.append(prob4[i])

auc4 = roc_auc_score(flags4, prob4)
    
fprx4, tprx4, _ = roc_curve(flags4.ravel(), prob4)

flags5 = np.load('flag_spanet_pretraining.npy')
prob5 = np.load('prob_spanet_pretraining.npy')

presig = []
prebkg = []
for i in range(len(flags5)):
    if flags5[i]==0:
        prebkg.append(prob5[i])
    else:
        presig.append(prob5[i])


auc5 = roc_auc_score(flags5, prob5)
    
fprx5, tprx5, _ = roc_curve(flags5.ravel(), prob5)

#flags6= np.load('flag_span_w_sig_tr_te_nopro.npy')
#prob6 = np.load('prob_span_w_sig_tr_te_nopro.npy')


#auc6 = roc_auc_score(flags6, prob6)
    
#fprx6, tprx6, _ = roc_curve(flags6.ravel(), prob6)

plt.plot(tprx4, 1-fprx4, linestyle = '-',label=r'{\textsc{SPANet}} Fine-tuning (AUC = '+str(round(auc4,3))+")",color = 'blue')
plt.plot(tprx3, 1-fprx3, linestyle = '--',label=r'{\textsc{SPANet}} BDT (AUC = '+str(round(auc3,3))+")",color = 'pink')
plt.plot(tprx5, 1-fprx5, linestyle = '-.',label=r'{\textsc{SPANet}} Pretraining (AUC = '+str(round(auc5,3))+")",color = 'green')
plt.plot(tprx2, 1-fprx2, linestyle = ':',label=r'PDNN BDT (AUC = '+str(round(auc2,3))+")",color = 'red')
#plt.plot(tprx6, 1-fprx6, linestyle = (0,(4, 1, 1, 1, 1, 1)),label='SPANet BDT no probabilities, AUC = '+str(round(auc6,3)),color = 'purple')
plt.plot(tprx, 1-fprx, linestyle = (0,(3, 1, 1, 1, 1, 1)) ,label=r'KLFitter BDT (AUC = '+str(round(auc1,3))+")",color = 'orange')
plt.legend(prop = {'size':25})
#plt.savefig('roc_combine.png')
plt.savefig('roc_combine.pdf')

# Save the model
norm = False
#nbin5
nbins = 5
klfbin5 = np.array([0, 0.35, 0.47, 0.56, 0.64, 1])
dnnbin5 = np.array([0, 0.35, 0.47, 0.56, 0.64, 1])
#spanbin5 = np.array([0, 0.210, 0.410, 0.61, 0.71, 1])、
spanbin5 = np.array([0, 0.26, 0.46, 0.59, 0.7, 1])
prebin5 = np.array([0, 0.45, 0.67, 0.78, 0.85, 1])
ftbin5 = np.array([0, 0.38, 0.61, 0.74, 0.84, 1])

plt.figure()
nklf, bins,patches = plt.hist(klfbkg,alpha=0.5,label='ttbar KLFitter',bins=klfbin5, density=norm)
mklf, bins, _ = plt.hist(klfsig, alpha=0.5, label="ttH KLFitter",bins=klfbin5,density=norm)
ntotklf=sum(nklf)
mtotklf=sum(mklf)
np.save('klf_sig_bin_nbins'+str(nbins),mklf)
np.save('klf_bkg_bin_nbins'+str(nbins),nklf)

ndnn, bins,patches = plt.hist(dnnbkg,alpha=0.5,label='ttbar PDNN',bins=dnnbin5, density=norm)
mdnn, bins, _ = plt.hist(dnnsig, alpha=0.5, label="ttH PDNN",bins=dnnbin5,density=norm)
ntotdnn=sum(ndnn)
mtotdnn=sum(mdnn)
np.save('dnn_sig_bin_nbins'+str(nbins),mdnn)
np.save('dnn_bkg_bin_nbins'+str(nbins),ndnn)

nspan, bins,patches = plt.hist(spanbkg,alpha=0.5,label='ttbar SPANet',bins=spanbin5, density=norm)
mspan, bins, _ = plt.hist(spansig, alpha=0.5, label="ttH SPANet",bins=spanbin5,density=norm)
ntotspan=sum(nspan)
mtotspan=sum(mspan)
np.save('span_sig_bin_nbins'+str(nbins),mspan)
np.save('span_bkg_bin_nbins'+str(nbins),nspan)

npre, bins,patches = plt.hist(prebkg,alpha=0.5,label='ttbar Pretraining',bins=prebin5, density=norm)
mpre, bins, _ = plt.hist(presig, alpha=0.5, label="ttH Pretraining",bins=prebin5,density=norm)
ntotpre=sum(npre)
mtotpre=sum(mpre)
np.save('pre_sig_bin_nbins'+str(nbins),mpre)
np.save('pre_bkg_bin_nbins'+str(nbins),npre)

nft, bins,patches = plt.hist(ftbkg,alpha=0.5,label='ttbar Fine-tuning',bins=ftbin5, density=norm)
mft, bins, _ = plt.hist(ftsig, alpha=0.5, label="ttH Fine-tuning",bins=ftbin5,density=norm)
ntotft=sum(nft)
mtotft=sum(mft)
np.save('ft_sig_bin_nbins'+str(nbins),mft)
np.save('ft_bkg_bin_nbins'+str(nbins),nft)

klfbincentre5 = getbincentre(klfbin5)
dnnbincentre5 = getbincentre(dnnbin5)
spanbincentre5 = getbincentre(spanbin5)
prebincentre5 = getbincentre(prebin5)
ftbincentre5 = getbincentre(ftbin5)

plt.figure(figsize=(12, 10))
metric=wasserstein_distance(klfsig, klfbkg)
plt.hist(klfbincentre5, bins=klfbin5, weights=nklf/ntotklf,alpha=0.5, label='ttbar KLFitter BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'orange')
plt.hist(klfbincentre5, bins=klfbin5, weights=mklf/mtotklf,alpha=0.5, label='ttH KLFitter BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'orange')
metric=wasserstein_distance(dnnsig, dnnbkg)
plt.hist(dnnbincentre5, bins=dnnbin5, weights=ndnn/ntotdnn,alpha=0.5, label='ttbar PDNN BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'red')
plt.hist(dnnbincentre5, bins=dnnbin5, weights=mdnn/mtotdnn,alpha=0.5, label='ttH PDNN BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'red')
metric=wasserstein_distance(spansig, spanbkg)
plt.hist(spanbincentre5, bins=spanbin5, weights=nspan/ntotspan,alpha=0.5, label='ttbar SPANet BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'pink')
plt.hist(spanbincentre5, bins=spanbin5, weights=mspan/mtotspan,alpha=0.5, label='ttH SPANet BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'pink')
metric=wasserstein_distance(presig, prebkg)
plt.hist(prebincentre5, bins=prebin5, weights=npre/ntotpre,alpha=0.5, label='ttbar Pretraining BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'green')
plt.hist(prebincentre5, bins=prebin5, weights=mpre/mtotpre,alpha=0.5, label='ttH Pretraining BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'green')
metric=wasserstein_distance(ftsig, ftbkg)
plt.hist(ftbincentre5, bins=ftbin5, weights=nft/ntotft,alpha=0.5, label='ttbar Fine-tuning BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'blue')
plt.hist(ftbincentre5, bins=ftbin5, weights=mft/mtotft,alpha=0.5, label='ttH Fine-tuning BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'blue')

plt.errorbar(klfbincentre5, nklf/ntotklf, np.sqrt(nklf)/ntotklf, fmt='none', ecolor='orange')
plt.errorbar(klfbincentre5, nklf/ntotklf, np.sqrt(mklf)/mtotklf, fmt='none', ecolor='orange')
plt.errorbar(dnnbincentre5, ndnn/ntotdnn, np.sqrt(ndnn)/ntotdnn, fmt='none', ecolor='red')
plt.errorbar(dnnbincentre5, ndnn/ntotdnn, np.sqrt(mdnn)/mtotdnn, fmt='none', ecolor='red')
plt.errorbar(spanbincentre5, nspan/ntotspan, np.sqrt(nspan)/ntotspan, fmt='none', ecolor='pink')
plt.errorbar(spanbincentre5, nspan/ntotspan, np.sqrt(mspan)/mtotspan, fmt='none', ecolor='pink')
plt.errorbar(prebincentre5, npre/ntotpre, np.sqrt(npre)/ntotpre, fmt='none', ecolor='green')
plt.errorbar(prebincentre5, npre/ntotpre, np.sqrt(mpre)/mtotpre, fmt='none', ecolor='green')
plt.errorbar(ftbincentre5, nft/ntotft, np.sqrt(nft)/ntotft, fmt='none', ecolor='blue')
plt.errorbar(ftbincentre5, nft/ntotft, np.sqrt(mft)/mtotft, fmt='none', ecolor='blue')


plt.legend(loc='upper left', bbox_to_anchor=(0.05, 0.95))
plt.savefig("classifieroutput_normalized"+str(nbins)+".png")
#nbin10
nbins = 10
klfbin10 = np.array([0, 0.15, 0.29, 0.38, 0.45, 0.5, 0.55, 0.59, 0.63, 0.68, 1])
dnnbin10 = np.array([0, 0.15, 0.29, 0.38, 0.45, 0.5, 0.55, 0.59, 0.64, 0.68, 1])
spanbin10 = np.array([0, 0.08, 0.22, 0.34, 0.44, 0.51, 0.58, 0.64, 0.7, 0.77, 1])
prebin10 = np.array([0, 0.04, 0.32, 0.5, 0.62, 0.7, 0.76, 0.8, 0.84, 0.88, 1])
ftbin10 = np.array([0, 0.14, 0.35, 0.5, 0.6, 0.67, 0.73, 0.79, 0.84, 0.89, 1])

klfbincentre10 = getbincentre(klfbin10)
dnnbincentre10 = getbincentre(dnnbin10)
spanbincentre10 = getbincentre(spanbin10)
prebincentre10 = getbincentre(prebin10)
ftbincentre10 = getbincentre(ftbin10)

plt.figure()
nklf, bins,patches = plt.hist(klfbkg,alpha=0.5,label='ttbar KLFitter',bins=klfbin10, density=norm)
mklf, bins, _ = plt.hist(klfsig, alpha=0.5, label="ttH KLFitter",bins=klfbin10,density=norm)
ntotklf=sum(nklf)
mtotklf=sum(mklf)
np.save('klf_sig_bin_nbins'+str(nbins),mklf)
np.save('klf_bkg_bin_nbins'+str(nbins),nklf)

ndnn, bins,patches = plt.hist(dnnbkg,alpha=0.5,label='ttbar PDNN',bins=dnnbin10, density=norm)
mdnn, bins, _ = plt.hist(dnnsig, alpha=0.5, label="ttH PDNN",bins=dnnbin10,density=norm)
ntotdnn=sum(ndnn)
mtotdnn=sum(mdnn)
np.save('dnn_sig_bin_nbins'+str(nbins),mdnn)
np.save('dnn_bkg_bin_nbins'+str(nbins),ndnn)

nspan, bins,patches = plt.hist(spanbkg,alpha=0.5,label='ttbar SPANet',bins=spanbin10, density=norm)
mspan, bins, _ = plt.hist(spansig, alpha=0.5, label="ttH SPANet",bins=spanbin10,density=norm)
ntotspan=sum(nspan)
mtotspan=sum(mspan)
np.save('span_sig_bin_nbins'+str(nbins),mspan)
np.save('span_bkg_bin_nbins'+str(nbins),nspan)

npre, bins,patches = plt.hist(prebkg,alpha=0.5,label='ttbar Pretraining',bins=prebin10, density=norm)
mpre, bins, _ = plt.hist(presig, alpha=0.5, label="ttH Pretraining",bins=prebin10,density=norm)
ntotpre=sum(npre)
mtotpre=sum(mpre)
np.save('pre_sig_bin_nbins'+str(nbins),mpre)
np.save('pre_bkg_bin_nbins'+str(nbins),npre)

nft, bins,patches = plt.hist(ftbkg,alpha=0.5,label='ttbar Fine-tuning',bins=ftbin10, density=norm)
mft, bins, _ = plt.hist(ftsig, alpha=0.5, label="ttH Fint-tuning",bins=ftbin10,density=norm)
ntotft=sum(nft)
mtotft=sum(mft)
np.save('ft_sig_bin_nbins'+str(nbins),mft)
np.save('ft_bkg_bin_nbins'+str(nbins),nft)

klfbincentre10 = getbincentre(klfbin10)
dnnbincentre10 = getbincentre(dnnbin10)
spanbincentre10 = getbincentre(spanbin10)
prebincentre10 = getbincentre(prebin10)
ftbincentre10 = getbincentre(ftbin10)

plt.figure(figsize=(12, 10))
metric=wasserstein_distance(klfsig, klfbkg)
plt.hist(klfbincentre10, bins=klfbin10, weights=nklf/ntotklf,alpha=0.5, label='ttbar KLFitter BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'orange')
plt.hist(klfbincentre10, bins=klfbin10, weights=mklf/mtotklf,alpha=0.5, label='ttH KLFitter BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'orange')
metric=wasserstein_distance(dnnsig, dnnbkg)
plt.hist(dnnbincentre10, bins=dnnbin10, weights=ndnn/ntotdnn,alpha=0.5, label='ttbar PDNN BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'red')
plt.hist(dnnbincentre10, bins=dnnbin10, weights=mdnn/mtotdnn,alpha=0.5, label='ttH PDNN BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'red')
metric=wasserstein_distance(spansig, spanbkg)
plt.hist(spanbincentre10, bins=spanbin10, weights=nspan/ntotspan,alpha=0.5, label='ttbar SPANet BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'pink')
plt.hist(spanbincentre10, bins=spanbin10, weights=mspan/mtotspan,alpha=0.5, label='ttH SPANet BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'pink')
metric=wasserstein_distance(presig, prebkg)
plt.hist(prebincentre10, bins=prebin10, weights=npre/ntotpre,alpha=0.5, label='ttbar Pretraining BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'green')
plt.hist(prebincentre10, bins=prebin10, weights=mpre/mtotpre,alpha=0.5, label='ttH Pretraining BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'green')
metric=wasserstein_distance(ftsig, ftbkg)
plt.hist(ftbincentre10, bins=ftbin10, weights=nft/ntotft,alpha=0.5, label='ttbar Fine-tuning BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'blue')
plt.hist(ftbincentre10, bins=ftbin10, weights=mft/mtotft,alpha=0.5, label='ttH Fine-tuning BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'blue')

plt.errorbar(klfbincentre10, nklf/ntotklf, np.sqrt(nklf)/ntotklf, fmt='none', ecolor='orange')
plt.errorbar(klfbincentre10, nklf/ntotklf, np.sqrt(mklf)/mtotklf, fmt='none', ecolor='orange')
plt.errorbar(dnnbincentre10, ndnn/ntotdnn, np.sqrt(ndnn)/ntotdnn, fmt='none', ecolor='red')
plt.errorbar(dnnbincentre10, ndnn/ntotdnn, np.sqrt(mdnn)/mtotdnn, fmt='none', ecolor='red')
plt.errorbar(spanbincentre10, nspan/ntotspan, np.sqrt(nspan)/ntotspan, fmt='none', ecolor='pink')
plt.errorbar(spanbincentre10, nspan/ntotspan, np.sqrt(mspan)/mtotspan, fmt='none', ecolor='pink')
plt.errorbar(prebincentre10, npre/ntotpre, np.sqrt(npre)/ntotpre, fmt='none', ecolor='green')
plt.errorbar(prebincentre10, npre/ntotpre, np.sqrt(mpre)/mtotpre, fmt='none', ecolor='green')
plt.errorbar(ftbincentre10, nft/ntotft, np.sqrt(nft)/ntotft, fmt='none', ecolor='blue')
plt.errorbar(ftbincentre10, nft/ntotft, np.sqrt(mft)/mtotft, fmt='none', ecolor='blue')


plt.legend(loc='upper left', bbox_to_anchor=(0.05, 0.95))
plt.savefig("classifieroutput_normalized"+str(nbins)+".png")

nbins = 8
klfbin8 = np.array([0, 0.2, 0.34, 0.43, 0.5, 0.56, 0.61, 0.67, 1])
dnnbin8 = np.array([0, 0.2, 0.34, 0.43, 0.5, 0.56, 0.61, 0.67, 1])
spanbin8 = np.array([0, 0.11, 0.28, 0.42, 0.51, 0.59, 0.66, 0.74, 1])
prebin8 = np.array([0, 0.23, 0.48, 0.63, 0.72, 0.78, 0.83, 0.87, 1])
ftbin8 = np.array([0, 0.24, 0.46, 0.59, 0.68, 0.75, 0.82, 0.88, 1])

klfbincentre8 = getbincentre(klfbin8)
dnnbincentre8 = getbincentre(dnnbin8)
spanbincentre8 = getbincentre(spanbin8)
prebincentre8 = getbincentre(prebin8)
ftbincentre8 = getbincentre(ftbin8)

plt.figure()
nklf, bins,patches = plt.hist(klfbkg,alpha=0.5,label='ttbar KLFitter',bins=klfbin8, density=norm)
mklf, bins, _ = plt.hist(klfsig, alpha=0.5, label="ttH KLFitter",bins=klfbin8,density=norm)
ntotklf=sum(nklf)
mtotklf=sum(mklf)
np.save('klf_sig_bin_nbins'+str(nbins),mklf)
np.save('klf_bkg_bin_nbins'+str(nbins),nklf)

ndnn, bins,patches = plt.hist(dnnbkg,alpha=0.5,label='ttbar PDNN',bins=dnnbin8, density=norm)
mdnn, bins, _ = plt.hist(dnnsig, alpha=0.5, label="ttH PDNN",bins=dnnbin8,density=norm)
ntotdnn=sum(ndnn)
mtotdnn=sum(mdnn)
np.save('dnn_sig_bin_nbins'+str(nbins),mdnn)
np.save('dnn_bkg_bin_nbins'+str(nbins),ndnn)

nspan, bins,patches = plt.hist(spanbkg,alpha=0.5,label='ttbar SPANet',bins=spanbin8, density=norm)
mspan, bins, _ = plt.hist(spansig, alpha=0.5, label="ttH SPANet",bins=spanbin8,density=norm)
ntotspan=sum(nspan)
mtotspan=sum(mspan)
np.save('span_sig_bin_nbins'+str(nbins),mspan)
np.save('span_bkg_bin_nbins'+str(nbins),nspan)

npre, bins,patches = plt.hist(prebkg,alpha=0.5,label='ttbar Pretraining',bins=prebin8, density=norm)
mpre, bins, _ = plt.hist(presig, alpha=0.5, label="ttH Pretraining",bins=prebin8,density=norm)
ntotpre=sum(npre)
mtotpre=sum(mpre)
np.save('pre_sig_bin_nbins'+str(nbins),mpre)
np.save('pre_bkg_bin_nbins'+str(nbins),npre)

nft, bins,patches = plt.hist(ftbkg,alpha=0.5,label='ttbar Fine-tuning',bins=ftbin8, density=norm)
mft, bins, _ = plt.hist(ftsig, alpha=0.5, label="ttH Fint-tuning",bins=ftbin8,density=norm)
ntotft=sum(nft)
mtotft=sum(mft)
np.save('ft_sig_bin_nbins'+str(nbins),mft)
np.save('ft_bkg_bin_nbins'+str(nbins),nft)

klfbincentre8 = getbincentre(klfbin8)
dnnbincentre8 = getbincentre(dnnbin8)
spanbincentre8 = getbincentre(spanbin8)
prebincentre8 = getbincentre(prebin8)
ftbincentre8 = getbincentre(ftbin8)

plt.figure(figsize=(12, 8))
metric=wasserstein_distance(klfsig, klfbkg)
plt.hist(klfbincentre8, bins=klfbin8, weights=nklf/ntotklf,alpha=0.5, label='ttbar KLFitter BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'orange')
plt.hist(klfbincentre8, bins=klfbin8, weights=mklf/mtotklf,alpha=0.5, label='ttH KLFitter BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'orange')
metric=wasserstein_distance(dnnsig, dnnbkg)
plt.hist(dnnbincentre8, bins=dnnbin8, weights=ndnn/ntotdnn,alpha=0.5, label='ttbar PDNN BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'red')
plt.hist(dnnbincentre8, bins=dnnbin8, weights=mdnn/mtotdnn,alpha=0.5, label='ttH PDNN BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'red')
metric=wasserstein_distance(spansig, spanbkg)
plt.hist(spanbincentre8, bins=spanbin8, weights=nspan/ntotspan,alpha=0.5, label='ttbar SPANet BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'pink')
plt.hist(spanbincentre8, bins=spanbin8, weights=mspan/mtotspan,alpha=0.5, label='ttH SPANet BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'pink')
metric=wasserstein_distance(presig, prebkg)
plt.hist(prebincentre8, bins=prebin8, weights=npre/ntotpre,alpha=0.5, label='ttbar Pretraining BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'green')
plt.hist(prebincentre8, bins=prebin8, weights=mpre/mtotpre,alpha=0.5, label='ttH Pretraining BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'green')
metric=wasserstein_distance(ftsig, ftbkg)
plt.hist(ftbincentre8, bins=ftbin8, weights=nft/ntotft,alpha=0.5, label='ttbar Fine-tuning BDT'+str(nbins)+'bins',histtype="step",linestyle = '--',color = 'blue')
plt.hist(ftbincentre8, bins=ftbin8, weights=mft/mtotft,alpha=0.5, label='ttH Fine-tuning BDT'+str(nbins)+'bins, EMD='+str(round(metric,3)),histtype="step",color = 'blue')

plt.errorbar(klfbincentre8, nklf/ntotklf, np.sqrt(nklf)/ntotklf, fmt='none', ecolor='orange')
plt.errorbar(klfbincentre8, nklf/ntotklf, np.sqrt(mklf)/mtotklf, fmt='none', ecolor='orange')
plt.errorbar(dnnbincentre8, ndnn/ntotdnn, np.sqrt(ndnn)/ntotdnn, fmt='none', ecolor='red')
plt.errorbar(dnnbincentre8, ndnn/ntotdnn, np.sqrt(mdnn)/mtotdnn, fmt='none', ecolor='red')
plt.errorbar(spanbincentre8, nspan/ntotspan, np.sqrt(nspan)/ntotspan, fmt='none', ecolor='pink')
plt.errorbar(spanbincentre8, nspan/ntotspan, np.sqrt(mspan)/mtotspan, fmt='none', ecolor='pink')
plt.errorbar(prebincentre8, npre/ntotpre, np.sqrt(npre)/ntotpre, fmt='none', ecolor='green')
plt.errorbar(prebincentre8, npre/ntotpre, np.sqrt(mpre)/mtotpre, fmt='none', ecolor='green')
plt.errorbar(ftbincentre8, nft/ntotft, np.sqrt(nft)/ntotft, fmt='none', ecolor='blue')
plt.errorbar(ftbincentre8, nft/ntotft, np.sqrt(mft)/mtotft, fmt='none', ecolor='blue')


plt.legend(loc='upper left', bbox_to_anchor=(0.05, 0.95))
plt.savefig("classifieroutput_normalized"+str(nbins)+".png")

