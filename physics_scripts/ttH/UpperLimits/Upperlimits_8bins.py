from pathlib import Path
from math import sqrt

import numpy as np
import matplotlib.pyplot as plt

import pyhf
import pyhf.readxml
from pyhf.contrib.viz import brazil

SMxec = 0.507
unit = 'pb'
nbins=8
kttbb = 1.5

def loadData(inputs):

    if inputs=="KLFitterBDT":
        bkgBinning = np.load('klf_bkg_bin_nbins'+str(nbins)+'.npy')
        sigBinning = np.load('klf_sig_bin_nbins'+str(nbins)+'.npy')

    if inputs=="DNNBDT":
        bkgBinning = np.load('dnn_bkg_bin_nbins'+str(nbins)+'.npy')
        sigBinning = np.load('dnn_sig_bin_nbins'+str(nbins)+'.npy')

    if inputs=="SPANetBDT":
        bkgBinning = np.load('span_bkg_bin_nbins'+str(nbins)+'.npy')
        sigBinning = np.load('span_sig_bin_nbins'+str(nbins)+'.npy')

    if inputs=="Pretraining":
        bkgBinning = np.load('pre_bkg_bin_nbins'+str(nbins)+'.npy')
        sigBinning = np.load('pre_sig_bin_nbins'+str(nbins)+'.npy')

    if inputs=="Finetunning":
       bkgBinning = np.load('ft_bkg_bin_nbins'+str(nbins)+'.npy')
       sigBinning = np.load('ft_sig_bin_nbins'+str(nbins)+'.npy')


    return bkgBinning, sigBinning

def gettingLimits(inputs):

    bkgBinning, sigBinning = loadData(inputs)
    sigcounting = 502  # LUM = 139 fb^-1
    bkgcounting = 4337
    #sigcounting = 213 + 113 + 59.9 + 13.9 + 3.09 + 35.1 + 8.5 #pre-fit
    #bkgcounting = 3160 + 1530 + 720 + 215 + 55 + 246 + 55
    #bkgUncercounting = sqrt(490*490 + 210*210 + 130*130 + 61*61 + 27*27 + 49*49 + 24*24)
    #sigcounting = 76 + 40 + 21 + 4.8 + 1 + 12 + 3 #post-fit
    #bkgcounting = 4530 + 2050 + 865 + 240 + 44.7 + 295 + 52.5
    #bkgUncercounting = sqrt(160*160 + 80*80 + 47*47 + 22*22 + 8.5*8.5 + 23*23 + 9.3*9.3)
    kttbb = 1.5
    if nbins==5:
        bkgUncerrate = 0.1
        #bkgUncerrate = bkgUncercounting/bkgcounting
        sigcounting = sigcounting*140/139 # LUM = 140 fb^-1 Run2 Style
        bkgcounting = bkgcounting*140/139
    else:
        bkgUncerrate = 0.07
        sigcounting = sigcounting*300/139 # LUM = 300 fb^-1 Run3 Style
        bkgcounting = bkgcounting*300/139
    
    bkgUncerBinning = []
    sigBinningNormalized = []
    bkgBinningNormalized  = []
    dataBinningNormalized = []
    alpha = 0.05
    for i in range(0,len(bkgBinning)):
       bkgUncerBinning.append(bkgUncerrate*kttbb*bkgBinning[i]*bkgcounting/sum(bkgBinning))
       sigBinningNormalized.append(sigBinning[i]*sigcounting/sum(sigBinning))
       bkgBinningNormalized.append(kttbb*bkgBinning[i]*bkgcounting/sum(bkgBinning)) 
       dataBinningNormalized.append(kttbb*bkgBinning[i]*bkgcounting/sum(bkgBinning)+sigBinning[i]*sigcounting/sum(sigBinning))

    

    pyhf.set_backend("numpy")
    model = pyhf.simplemodels.uncorrelated_background(
       signal=sigBinningNormalized, bkg=bkgBinningNormalized, bkg_uncertainty= bkgUncerBinning
    )
    data = dataBinningNormalized + model.config.auxdata

    poi_vals = np.linspace(0, 8, 800)


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
    fig.savefig(inputs+str(nbins)+".png")

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