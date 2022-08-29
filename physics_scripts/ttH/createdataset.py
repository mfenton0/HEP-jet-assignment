from fileinput import filename
from math import cos, fabs, sin
from select import select
from traceback import print_tb
from types import MethodType
import h5py
import numpy as np
import ROOT
import math

trainsize=173114
testsize=30000
nosig = 0
filename = 'tth_with_spanet_with_signal.h5'
if nosig == 1:
  trainsize=173108
  filename = 'tth_with_spanet_no_signal.h5'

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

    bhadindex = file["spanet"]["right_target"]["b"][:trainsize]
    q1index = file["spanet"]["right_target"]["q1"][:trainsize]
    q2index = file["spanet"]["right_target"]["q2"][:trainsize]
    blepindex = file["spanet"]["left_target"]["b"][:trainsize]
    b1index = file["spanet"]["higgs_target"]["b1"][:trainsize]
    b2index = file["spanet"]["higgs_target"]["b2"][:trainsize]

    flags = file["spanet"]["signal"]["target"][:trainsize]
    spanet_probability = file["spanet"]["signal"]["probability"][:trainsize]
    assignment_probability = file["spanet"]["higgs_target"]["assignment_probability"][:trainsize]
    detection_probability = file["spanet"]["higgs_target"]["detection_probability"][:trainsize]
    marginal_probability = file["spanet"]["higgs_target"]["marginal_probability"][:trainsize]


    Ravgbb=np.zeros(trainsize)
    Rmaxptbb=np.zeros(trainsize)
    massminRbb=np.zeros(trainsize)
    Nbbhiggs30=np.zeros(trainsize)
    Rminlbb=np.zeros(trainsize)
    massminRjj=np.zeros(trainsize)
    etamaxjj=np.zeros(trainsize)
    HThad=np.zeros(trainsize)

    higgsm=np.zeros(trainsize)
    mhblep=np.zeros(trainsize)
    Rhbb=np.zeros(trainsize)
    Rhtt=np.zeros(trainsize)
    Rhtlep=np.zeros(trainsize)
    Rhbhad=np.zeros(trainsize)

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
        
        b1 = ROOT.TLorentzVector()
        b2 = ROOT.TLorentzVector()
        q1 = ROOT.TLorentzVector()
        q2 = ROOT.TLorentzVector()

        bhad = ROOT.TLorentzVector()
        blep = ROOT.TLorentzVector()

        Whad = ROOT.TLorentzVector()
        Wlep = ROOT.TLorentzVector()

        neutrino = ROOT.TLorentzVector()

        Thad = ROOT.TLorentzVector()
        Tlep = ROOT.TLorentzVector()

        higgs = ROOT.TLorentzVector()

        b1.SetPtEtaPhiM(pt[i][b1index[i]],eta[i][b1index[i]],phi[i][b1index[i]],mass[i][b1index[i]])
        b2.SetPtEtaPhiM(pt[i][b2index[i]],eta[i][b2index[i]],phi[i][b2index[i]],mass[i][b2index[i]])
        bhad.SetPtEtaPhiM(pt[i][bhadindex[i]],eta[i][bhadindex[i]],phi[i][bhadindex[i]],mass[i][bhadindex[i]])
        blep.SetPtEtaPhiM(pt[i][blepindex[i]],eta[i][blepindex[i]],phi[i][blepindex[i]],mass[i][blepindex[i]])
        q1.SetPtEtaPhiM(pt[i][q1index[i]],eta[i][q1index[i]],phi[i][q1index[i]],mass[i][q1index[i]])
        q2.SetPtEtaPhiM(pt[i][q2index[i]],eta[i][q2index[i]],phi[i][q2index[i]],mass[i][q2index[i]])

        higgs = b1 + b2
        Whad = q1 + q2
        Thad = Whad + bhad

        hhh = (80.385)*(80.385) - lepton.M()*lepton.M() + 2*metpt[i]*(math.cos(metphi[i])*lepton.Px()+math.sin(metphi[i])*lepton.Py())
        aaa = 4*(lepton.Pz()*lepton.Pz()) - 4*lepton.E()*lepton.E()
        bbb = 4*lepton.Pz()*hhh
        ccc = hhh*hhh - 4*metpt[i]*metpt[i]*lepton.E()*lepton.E()
        if (bbb*bbb - 4*aaa*ccc)>0 :
            solution1 = (math.sqrt(bbb*bbb - 4*aaa*ccc)-bbb)/(2*aaa)
            solution2 = (-math.sqrt(bbb*bbb - 4*aaa*ccc)-bbb)/(2*aaa)
            if fabs(solution1)<fabs(solution2):
                mez = solution1
            else:
                mez = solution2
        else:
            mez = -bbb/(2*aaa)

        mex=metpt[i]*math.cos(metphi[i])
        mey=metpt[i]*math.sin(metphi[i])

        neutrino.SetPxPyPzE(mex,mey,mez,math.sqrt(mex*mex+mey*mey+mez*mez))

        Wlep=neutrino + lepton
        Tlep=Wlep + blep

        higgsm[i] = higgs.M()
        hblep = higgs + blep
        mhblep[i] = hblep.M()
        Rhbb[i] = b1.DeltaR(b2,False)
        ttbar = Tlep + Thad
        Rhtt[i] = higgs.DeltaR(ttbar,False)
        Rhtlep[i] = higgs.DeltaR(Tlep,False)
        Rhbhad[i] = higgs.DeltaR(bhad,False)

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
        ("flags", "f4"),
        ("spanet_probability", "f4"),
        ("assignment_probability", "f4"),
        ("detection_probability", "f4"),
        ("marginal_probability", "f4"),
        ("eventnumber", "f4"),
    ]
)

spanet_data = np.zeros(28780, dtype=DATA_STRUCTURE)
sel = 0
for fill in range(trainsize):
    if (bjetnum[fill]>0) & (jetcounts[fill]==8) :
      spanet_data[sel]["higgsm"] = higgsm[fill]
      spanet_data[sel]["Ravgbb"] = Ravgbb[fill]
      spanet_data[sel]["Rmaxptbb"] = Rmaxptbb[fill]
      spanet_data[sel]["etamaxjj"] = etamaxjj[fill]
      spanet_data[sel]["massminRbb"] = massminRbb[fill]
      spanet_data[sel]["massminRjj"] = massminRjj[fill]
      spanet_data[sel]["Nbbhiggs30"] = Nbbhiggs30[fill]
      spanet_data[sel]["HThad"] = HThad[fill]
      spanet_data[sel]["Rminlbb"] = Rminlbb[fill]
      spanet_data[sel]["mhblep"] = mhblep[fill]
      spanet_data[sel]["Rhbb"] = Rhbb[fill]
      spanet_data[sel]["Rhtt"] = Rhtt[fill]
      spanet_data[sel]["Rhtlep"] = Rhtlep[fill]
      spanet_data[sel]["Rhbhad"] = Rhbhad[fill]
      spanet_data[sel]["flags"] = flags[fill]
      spanet_data[sel]["spanet_probability"] = spanet_probability[fill]
      spanet_data[sel]["assignment_probability"] = assignment_probability[fill]
      spanet_data[sel]["detection_probability"] = detection_probability[fill]
      spanet_data[sel]["marginal_probability"] = marginal_probability[fill]
      spanet_data[sel]["eventnumber"] = fill
      
      sel = sel + 1
      #if flags[fill]==1:
      #  if b1index[fill]==b2index[fill] :
      #      print("Error in prediction") 
    

with h5py.File("spanet_features.h5", "w") as file:
    file.create_dataset("events", data=spanet_data, dtype=DATA_STRUCTURE)







    

print(sel)