from fileinput import filename
from math import cos, fabs, sin
from select import select
from traceback import print_tb
from types import MethodType
import h5py
import numpy as np
import ROOT
import math
from glob import glob
      
DATA_STRUCTURE = np.dtype(
      [
        ("Ravgbb", "f4"),
        ("Rmaxptbb", "f4"),
        ("etamaxjj", "f4"),
        ("massminRbb", "f4"),
        ("massminRjj", "f4"),
        ("Nbbhiggs30", "f4"),
        ("HThad", "f4"),
        ("Rminlbb", "f4"),
      ]
  )


def updatefeature(filename : str):
   with h5py.File(filename, 'r') as fileA:
    # Load in the raw data
      sizesearch = fileA["jet_features"]["pt"][:].astype(np.float32)
   trainsize = len(sizesearch)
   print(filename)
   with h5py.File(filename, 'r') as file:
    # Load in the raw data
      pt = file["jet_features"]["pt"][:trainsize].astype(np.float32)
      eta = file["jet_features"]["eta"][:trainsize].astype(np.float32)
      phi = file["jet_features"]["phi"][:trainsize].astype(np.float32)
      mass = file["jet_features"]["mass"][:trainsize].astype(np.float32)
      btag = file["jet_features"]["btag"][:trainsize].astype(np.float32)

      leptonpt = file["lepton_features"]["pt"][:trainsize].astype(np.float32)
      leptoneta = file["lepton_features"]["eta"][:trainsize].astype(np.float32)
      leptonphi = file["lepton_features"]["phi"][:trainsize].astype(np.float32)
      leptonmass = file["lepton_features"]["mass"][:trainsize].astype(np.float32)

      metpt = file["met_features"]["MET"][:trainsize].astype(np.float32)
      metphi = file["met_features"]["phi"][:trainsize].astype(np.float32)

      Ravgbb=np.zeros(trainsize)
      Rmaxptbb=np.zeros(trainsize)
      massminRbb=np.zeros(trainsize)
      Nbbhiggs30=np.zeros(trainsize)
      Rminlbb=np.zeros(trainsize)
      massminRjj=np.zeros(trainsize)
      etamaxjj=np.zeros(trainsize)
      HThad=np.zeros(trainsize)


      bjetnum=np.zeros(trainsize)
      jetcounts=[]
      for i in range(trainsize):
          jets=[]
          bindices=[]
          lepton=ROOT.TLorentzVector()
          lepton.SetPtEtaPhiM(leptonpt[i],leptoneta[i],leptonphi[i],leptonmass[i])
          numofjet = 0
          for j in range(len(pt[i])):
              if pt[i][j]>0 :
                  jet1=ROOT.TLorentzVector()
                  jet1.SetPtEtaPhiM(pt[i][j],eta[i][j],phi[i][j],mass[i][j])
                  jets.append(jet1)
                  numofjet=numofjet+1
                  if btag[i][j]>0:
                      bindices.append(j)
          jetcounts.append(numofjet)
          bjetnum[i]=len(bindices)
          sumRbb=[]
          jetsumbb0=jets[bindices[0]] + jets[bindices[1]]
          jetsumbbptmax = jetsumbb0.Pt()
          Rmaxptbb[i]=jets[bindices[0]].DeltaR(jets[bindices[1]],False)
          deltaRbbmin=jets[bindices[0]].DeltaR(jets[bindices[1]],False)
          massminRbb[i]=jetsumbb0.M()
          Nbbhiggs30[i]=0
          Rminlbb[i]=jetsumbb0.DeltaR(lepton,False)
          for i1 in range(len(bindices)):
              for j1 in range(len(bindices)):
                  if j1>i1:
                      deltaRij=jets[bindices[i1]].DeltaR(jets[bindices[j1]],False)
                      jetsumbb=jets[bindices[i1]] + jets[bindices[j1]]
                      jetsumbbpt = jetsumbb.Pt()
                      jetsumbbmass = jetsumbb.M()
                      if jetsumbbpt>jetsumbbptmax :
                          jetsumbbptmax=jetsumbbpt
                          Rmaxptbb[i]=deltaRij
                      if deltaRij<deltaRbbmin:
                          deltaRbbmin=deltaRij
                          massminRbb[i]=jetsumbb.M()
                          Rminlbb[i]=jetsumbb.DeltaR(lepton,False)
                      if (jetsumbbmass>95) & (jetsumbbmass<155):
                          Nbbhiggs30[i]=Nbbhiggs30[i]+1
                      sumRbb.append(deltaRij)
        
          Ravgbb[i]=sum(sumRbb)/len(sumRbb)
         
          jetsumjj0=jets[1] + jets[2]
          deltaRjjmin=jets[1].DeltaR(jets[2],False)
          massminRjj[i]=jetsumjj0.M()
          etamaxjj[i]=fabs(jets[1].Eta()-jets[2].Eta())
          for i2 in range(numofjet):
              for j2 in range(numofjet):
                  if j2>i2:
                      deltajjRij=jets[i2].DeltaR(jets[j2],False)
                      jetsumjj=jets[i2]+jets[j2]
                      if deltajjRij<deltaRjjmin:
                          deltaRjjmin=deltajjRij
                          massminRjj[i]=jetsumjj.M()
                      etajj=fabs(jets[i2].Eta()-jets[j2].Eta())
                      if etajj>etamaxjj[i]:
                          etamaxjj[i]=etajj

          jetpt=[]
          for m in range(numofjet):
              jetpt.append(jets[m].Pt())
          HThad[i]=sum(jetpt)
        
  

      spanet_data = np.zeros(trainsize, dtype=DATA_STRUCTURE)
      sel = 0
      for fill in range(trainsize):
        if (bjetnum[fill]>0) & (jetcounts[fill]>5) :
          spanet_data[sel]["Ravgbb"] = Ravgbb[fill]
          spanet_data[sel]["Rmaxptbb"] = Rmaxptbb[fill]
          spanet_data[sel]["etamaxjj"] = etamaxjj[fill]
          spanet_data[sel]["massminRbb"] = massminRbb[fill]
          spanet_data[sel]["massminRjj"] = massminRjj[fill]
          spanet_data[sel]["Nbbhiggs30"] = Nbbhiggs30[fill]
          spanet_data[sel]["HThad"] = HThad[fill]
          spanet_data[sel]["Rminlbb"] = Rminlbb[fill]
      
          sel = sel + 1

   return spanet_data

def main():
    input_folder = 'event_record_ttbb_lep_incl_CMS_jetR05_bPt20'
    all_files = list(glob(f"{input_folder}/*.h5"))
    for filename in all_files: 
      featureupdate = updatefeature(filename)
      with h5py.File(filename, "a") as file:
          file.create_dataset("KinematicVaribles", data=featureupdate, dtype=DATA_STRUCTURE)

if __name__ == "__main__":
    main()
