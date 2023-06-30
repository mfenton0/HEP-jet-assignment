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
kttbb = 1.5

def loadData(inputs):

    if inputs=="KLFitterBDT":
        bkgBinning = np.load('klf_bkg_bin_nbins'+str(nbins)+'.npy')
        sigBinning = np.load('klf_sig_bin_nbins'+str(nbins)+'.npy')

    if inputs=="DNNBDT":
        bkgBinning = np.load('dnn_bkg_bin_nbins'+str(nbins)+'.npy')
        sigBinning = np.load('dnn_sig_bin_nbins'+str(nbins)+'.npy')

    if inputs=="SPANetBDT":
        bkgBinning = np.load('span_bkg_bin_noprob_nbins'+str(nbins)+'.npy')
        sigBinning = np.load('span_sig_bin_noprob_nbins'+str(nbins)+'.npy')

    if inputs=="Pretraining":
        bkgBinning = np.load('span_bkg_bin_pretraining_nbins'+str(nbins)+'.npy')
        sigBinning = np.load('span_sig_bin_pretraining_nbins'+str(nbins)+'.npy')

    if inputs=="Finetunning":
       bkgBinning = np.load('span_bkg_bin_finetuning_nbin'+str(nbins)+'.npy')
       sigBinning = np.load('span_sig_bin_finetuning_nbin'+str(nbins)+'.npy')

    return bkgBinning, sigBinning

def gettingLimits(inputs):

    bkgBinning, sigBinning = loadData(inputs)
    sigcounting = 502  # LUM = 139 fb^-1
    bkgcounting = 4337
    if nbins==5:
        bkgUncerrate = 0.5
        sigcounting = sigcounting*140/139 # LUM = 140 fb^-1 Run2 Style
        bkgcounting = bkgcounting*140/139
    else:
        bkgrate = 0.3
        sigcounting = sigcounting*300/139 # LUM = 300 fb^-1 Run3 Style
        bkgcounting = bkgcounting*300/139
    
    bkgUncerBinning = []
    sigBinningNormalized = []
    bkgBinningNormalized  = []
    dataBinningNormalized = []
    alpha = 0.05
    for i in range(0,len(bkgBinning)):
       bkgUncerBinning.append(bkgUncerrate*sqrt(bkgBinning[i])*bkgcounting/sum(bkgBinning))
       sigBinningNormalized.append(sigBinning[i]*sigcounting/sum(sigBinning))
       bkgBinningNormalized.append(kttbb*bkgBinning[i]*bkgcounting/sum(bkgBinning)) 
       dataBinningNormalized.append(kttbb*bkgBinning[i]*bkgcounting/sum(bkgBinning)+sigBinning[i]*sigcounting/sum(sigBinning))

    pyhf.set_backend("numpy")
    model = pyhf.simplemodels.uncorrelated_background(
       signal=sigBinningNormalized, bkg=bkgBinningNormalized, bkg_uncertainty= bkgUncerBinning
    )
    data = dataBinningNormalized + model.config.auxdata

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
    plt.title(inputs+", expected limits:"+str(np.round(exp_limits[2],3))+" bkg uncertainty: 50% (300fb-1)")
    fig.savefig(inputs+".png")

    asymp_calc = pyhf.infer.calculators.AsymptoticCalculator(
       data, model, test_stat='qtilde'
    )
    teststat = asymp_calc.teststatistic(poi_test=obs_limit)
    print(inputs+":"+str(np.round(exp_limits[2],3))+" xsec:"+str(np.round((np.round(exp_limits[2],3)*SMxec),3))+"(pb) "+ f'qtilde = {teststat}')

def main():
    gettingLimits("KLFitterBDT")
    gettingLimits("DNNBDT")
    gettingLimits("SPANetBDT")
    gettingLimits("Pretraining")
    gettingLimits("Finetunning")

main()