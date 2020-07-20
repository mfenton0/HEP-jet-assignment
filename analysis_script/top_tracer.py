import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time, sys, os

PID_TOP = 6
PID_TOP_BAR = -6
STATUS = 22

class particle_properties():
    def __init__(self, dataset, pidx):
        self.pid = dataset.iloc[pidx, 6]
        self.pt = dataset.iloc[pidx, 7]
        self.eta = dataset.iloc[pidx, 8]
        self.phi = dataset.iloc[pidx, 9]
        self.mass = dataset.iloc[pidx, 10]

def shift_particle_tracing(dataset, PID_d, idx):
    if (dataset.iloc[idx,6] == PID_d):
        return dataset.iloc[idx,4]

def particle_tracing(dataset, PID):
    print("Start particle tracing!")
    for i in range(len(dataset)):
        if(dataset.iloc[i,1] == STATUS and dataset.iloc[i,6] == PID ): 
            daughter_index = dataset.iloc[i,0]
            #print(daughter_index)

    if( dataset.iloc[daughter_index,6] == PID ):
        shifted_particle_index = dataset.iloc[daughter_index, 4]
        #print(next_daughter_index)

    while dataset.iloc[shifted_particle_index,6] == PID:
            init_shifted_particle_index = shifted_particle_index
            shifted_particle_index = shift_particle_tracing(dataset, PID, init_shifted_particle_index)       
            #print("next_daughter_index: {0}".format(next_daughter_index))
    return init_shifted_particle_index

def deltaPhi(phi1,phi2):
    x = phi1-phi2
    while x>= np.pi: x -= np.pi*2.
    while x< -np.pi: x += np.pi*2.
    return x

def deltaR(eta1,phi1,eta2,phi2):
    return (deltaPhi(phi1,phi2)**2+(eta1-eta2)**2)**0.5
    
    
# def Match(jet, event, dR, _pid=1, _pT=2, _eta=3, _phi=4, _m=5, R = 0.5):
#     return (dR > deltaR(jet[1][0],jet[2][0],event[_eta],event[_phi]))

# def Match_dR(jet, event, _pid=5, _pT=6, _eta=7, _phi=8, _m=9, R = 0.5):
#     return deltaR(jet[1][0],jet[2][0],event[_eta],event[_phi])

def main():
    data = pd.read_csv("parsed_2020_07_20_11_20.txt")
    print("Length of data is:{0}, shape is:{1}".format(len(data), data.shape))
    
    patron_1_idx = particle_tracing(data, PID_TOP)
    patron_2_idx = particle_tracing(data, PID_TOP_BAR)
    print("Find partron before decay, index of patron 1 is {0}, index of patron 2 is {1}.".format(patron_1_idx, patron_2_idx))
    
    patron_1 = particle_properties(data, patron_1_idx)
    patron_2 = particle_properties(data, patron_2_idx)

    dR = deltaR(patron_1.eta, patron_1.phi, patron_2.eta, patron_2.phi)

    print("delta R of two patron is:{0:.3f}".format(dR))






if __name__ == '__main__':
    main()
