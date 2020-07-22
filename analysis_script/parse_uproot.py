#!/usr/bin/env python
# coding: utf-8



import uproot
import pandas as pd 
import numpy as np 
from particle_properties import particle_properties
from jet_properties import jet_properties



data  = uproot.open('./tag_1_delphes_events.root')['Delphes']
data.show()




particle = particle_properties(data)

Length = len(particle.event)
test_length = 10
#df0 = particle.dataframelize(0)
#df1 = particle.dataframelize(1)
#df2 = particle.dataframelize(2)




def shift_particle_tracing(dataset, PID_d, idx):
    if (dataset.iloc[idx,6] == PID_d):
        return dataset.iloc[idx,4]

def particle_tracing(dataset, PID, STATUS):

    for i in range(len(dataset)):
        if(dataset.iloc[i,1] == STATUS and dataset.iloc[i,6] == PID ): 
            daughter_index = int(dataset.iloc[i,0])
    if( dataset.iloc[daughter_index,6] == PID ):
        shifted_particle_index = dataset.iloc[daughter_index, 4]


    while dataset.iloc[shifted_particle_index,6] == PID:
            init_shifted_particle_index = shifted_particle_index
            shifted_particle_index = shift_particle_tracing(dataset, PID, init_shifted_particle_index)       

    dauthter_idx_1 = dataset.iloc[init_shifted_particle_index, 4]
    daughter_pid_1 = dataset.iloc[dauthter_idx_1, 6]

    dauthter_idx_2 = dataset.iloc[init_shifted_particle_index, 5]
    daughter_pid_2 = dataset.iloc[dauthter_idx_2, 6]

    return init_shifted_particle_index, dauthter_idx_1, daughter_pid_1, dauthter_idx_2, daughter_pid_2




PID_W_plus = 24 
PID_W_minus = -24
PID_DOWN = 1
PID_DOWN_VAR = -1
PID_UP = 2
PID_UP_BAR = -2
PID_STRANGE = 3
PID_STRANGE_BAR = -3
PID_CHARM = 4
PID_CHARM_BAR = -4
PID_BOTTOM = 5
PID_BOTTOM_BAR = -5
PID_TOP = 6
PID_TOP_BAR = -6

top_idx = np.zeros(len(particle.event))
top_daughter_idx_1 = np.zeros(len(particle.event))
top_daughter_pid_1 = np.zeros(len(particle.event))
top_daughter_idx_2 = np.zeros(len(particle.event))
top_daughter_pid_2 = np.zeros(len(particle.event))

top_bar_idx = np.zeros(len(particle.event))
top_bar_daughter_idx_1 = np.zeros(len(particle.event))
top_bar_daughter_pid_1 = np.zeros(len(particle.event))
top_bar_daughter_idx_2 = np.zeros(len(particle.event))
top_bar_daughter_pid_2 = np.zeros(len(particle.event))

quark_idx_1 = np.zeros(len(particle.event))

quark_idx_2 = np.zeros(len(particle.event))

quark_idx_3 = np.zeros(len(particle.event))

quark_idx_4 = np.zeros(len(particle.event))

quark_idx_5 = np.zeros(len(particle.event))

quark_idx_6 = np.zeros(len(particle.event))




for i in range(0,10):
    top_idx[i], top_daughter_idx_1[i], top_daughter_pid_1[i], top_daughter_idx_2[i], top_daughter_pid_2[i] = particle_tracing(particle.dataframelize(i), PID_TOP, 22)
    top_bar_idx[i], top_bar_daughter_idx_1[i], top_bar_daughter_pid_1[i], top_bar_daughter_idx_2[i], top_bar_daughter_pid_2[i] = particle_tracing(particle.dataframelize(i), PID_TOP_BAR, 22)




def quark_finder(dataset, mother_idx):

    def switcher(input_1, input_2, m_id):
        if m_id == 6:
            if input_1 == 24 :
                return  input_1, input_2
            else :
                return  input_2, input_1
        else :
            if input_1 == -24 :
                return  input_1, input_2
            else :
                return  input_2, input_1

    print("Mother index: {0}".format(mother_idx))

    mother_pid = dataset.iloc[int(mother_idx), 6]

    daughter_idx_1 = dataset.iloc[int(mother_idx), 4]
    daughter_idx_2 = dataset.iloc[int(mother_idx), 5]

    branch_daughter_W_idx, daught_b_idx = switcher(dataset.iloc[int(daughter_idx_1), 6], dataset.iloc[int(daughter_idx_2), 6],  mother_pid)

    status = dataset.iloc[int(branch_daughter_W_idx),1]
    print("Daughter's of t/t~'s status code: {0}".format(status))

    #find daughter of daughter 1
    while status != 23 :
        mother_idx = branch_daughter_W_idx
        daughter_idx_1_1 = dataset.iloc[int(mother_idx), 4]
        daughter_idx_1_2 = dataset.iloc[int(mother_idx), 5]
        status = dataset.iloc[int(daughter_idx_1), 1]
        print("2nd stage daughter status code: {0}".format(status))

    
    daughter_pid_1_1 = dataset.iloc[daughter_idx_1_1, 6]
    daughter_pid_1_2 = dataset.iloc[daughter_idx_1_2, 6]
    
    print("Daughter 1 index: {0}, daughter 2 index: {1}".format(daughter_idx_1_1, daughter_idx_1_2))
    print("Daughter 1 PID: {0}, daughter 2 PID: {1}".format(daughter_pid_1_1, daughter_pid_1_2))


    print( int(daught_b_idx), int(daughter_idx_1_1), int(daughter_idx_1_2))
    return int(daught_b_idx), int(daughter_idx_1_1), int(daughter_idx_1_2)




for i in range(0,10):
    quark_idx_1[i], quark_idx_2[i], quark_idx_3[i] = quark_finder(particle.dataframelize(i), top_daughter_idx_1[i])
    quark_idx_4[i], quark_idx_5[i], quark_idx_6[i] = quark_finder(particle.dataframelize(i), top_bar_daughter_idx_1[i])




quark_in_each_event = np.zeros([len(particle.event), 6, 6])
print(quark_in_each_event.shape)




for i in range(0,10):
    dataframe = particle.dataframelize(i)

    quark_in_each_event[i][0][0] = int(dataframe.iloc[int(quark_idx_1[i]),0])  #Index
    quark_in_each_event[i][0][1] = int(dataframe.iloc[int(quark_idx_1[i]),6])  #PID
    quark_in_each_event[i][0][2] = int(dataframe.iloc[int(quark_idx_1[i]),7])  #PT
    quark_in_each_event[i][0][3] = int(dataframe.iloc[int(quark_idx_1[i]),8])  #Eta
    quark_in_each_event[i][0][4] = int(dataframe.iloc[int(quark_idx_1[i]),9])  #Phi
    quark_in_each_event[i][0][5] = int(dataframe.iloc[int(quark_idx_1[i]),10])  #Mass

    quark_in_each_event[i][1][0] = int(dataframe.iloc[int(quark_idx_2[i]),0])  #Index
    quark_in_each_event[i][1][1] = int(dataframe.iloc[int(quark_idx_2[i]),6])  #PID
    quark_in_each_event[i][1][2] = int(dataframe.iloc[int(quark_idx_2[i]),7])  #PT
    quark_in_each_event[i][1][3] = int(dataframe.iloc[int(quark_idx_2[i]),8])  #Eta
    quark_in_each_event[i][1][4] = int(dataframe.iloc[int(quark_idx_2[i]),9])  #Phi
    quark_in_each_event[i][1][5] = int(dataframe.iloc[int(quark_idx_2[i]),10])  #Mass

    quark_in_each_event[i][2][0] = int(dataframe.iloc[int(quark_idx_3[i]),0])  #Index
    quark_in_each_event[i][2][1] = int(dataframe.iloc[int(quark_idx_3[i]),6])  #PID
    quark_in_each_event[i][2][2] = int(dataframe.iloc[int(quark_idx_3[i]),7])  #PT
    quark_in_each_event[i][2][3] = int(dataframe.iloc[int(quark_idx_3[i]),8])  #Eta
    quark_in_each_event[i][2][4] = int(dataframe.iloc[int(quark_idx_3[i]),9])  #Phi
    quark_in_each_event[i][2][5] = int(dataframe.iloc[int(quark_idx_3[i]),10])  #Mass

    quark_in_each_event[i][3][0] = int(dataframe.iloc[int(quark_idx_4[i]),0])  #Index
    quark_in_each_event[i][3][1] = int(dataframe.iloc[int(quark_idx_4[i]),6])  #PID
    quark_in_each_event[i][3][2] = int(dataframe.iloc[int(quark_idx_4[i]),7])  #PT
    quark_in_each_event[i][3][3] = int(dataframe.iloc[int(quark_idx_4[i]),8])  #Eta
    quark_in_each_event[i][3][4] = int(dataframe.iloc[int(quark_idx_4[i]),9])  #Phi
    quark_in_each_event[i][3][5] = int(dataframe.iloc[int(quark_idx_4[i]),10])  #Mass

    quark_in_each_event[i][4][0] = int(dataframe.iloc[int(quark_idx_5[i]),0])  #Index
    quark_in_each_event[i][4][1] = int(dataframe.iloc[int(quark_idx_5[i]),6])  #PID
    quark_in_each_event[i][4][2] = int(dataframe.iloc[int(quark_idx_5[i]),7])  #PT
    quark_in_each_event[i][4][3] = int(dataframe.iloc[int(quark_idx_5[i]),8])  #Eta
    quark_in_each_event[i][4][4] = int(dataframe.iloc[int(quark_idx_5[i]),9])  #Phi
    quark_in_each_event[i][4][5] = int(dataframe.iloc[int(quark_idx_5[i]),10])  #Mass

    quark_in_each_event[i][5][0] = int(dataframe.iloc[int(quark_idx_6[i]),0])  #Index
    quark_in_each_event[i][5][1] = int(dataframe.iloc[int(quark_idx_6[i]),6])  #PID
    quark_in_each_event[i][5][2] = int(dataframe.iloc[int(quark_idx_6[i]),7])  #PT
    quark_in_each_event[i][5][3] = int(dataframe.iloc[int(quark_idx_6[i]),8])  #Eta
    quark_in_each_event[i][5][4] = int(dataframe.iloc[int(quark_idx_6[i]),9])  #Phi
    quark_in_each_event[i][5][5] = int(dataframe.iloc[int(quark_idx_6[i]),10])  #Mass


jet = jet_properties(data)

print(len(jet.event), len(jet.pt), len(jet.pt[0]))




def deltaPhi(phi1,phi2):
    phi = phi1-phi2
    while phi >= np.pi: phi -= np.pi*2.
    while phi < -np.pi: phi += np.pi*2.
    return phi

def delta_R(eta1, phi1, eta2, phi2):
    return np.sqrt(deltaPhi(phi1,phi2)**2+(eta1-eta2)**2)




dR_patron_jet = np.zeros([len(particle.event), 4])
dR_patron_patron = np.zeros([len(particle.event), 6])




for i in range(0,10):
    dR_patron_patron[i][0] = delta_R(quark_in_each_event[i][0][3], quark_in_each_event[i][0][4], quark_in_each_event[i][1][3], quark_in_each_event[i][1][4] ) #dR(q1,q2)
    #dR_patron_patron[i][1] = delta_R(quark_in_each_event[i][0][3], quark_in_each_event[i][0][4], quark_in_each_event[i][2][3], quark_in_each_event[i][2][4] ) #dR(q1,q3)
    #dR_patron_patron[i][2] = delta_R(quark_in_each_event[i][0][3], quark_in_each_event[i][0][4], quark_in_each_event[i][3][3], quark_in_each_event[i][3][4] ) #dR(q1,q4)
    #dR_patron_patron[i][3] = delta_R(quark_in_each_event[i][1][3], quark_in_each_event[i][1][4], quark_in_each_event[i][2][3], quark_in_each_event[i][2][4] ) #dR(q2,q3)
    #dR_patron_patron[i][4] = delta_R(quark_in_each_event[i][1][3], quark_in_each_event[i][1][4], quark_in_each_event[i][3][3], quark_in_each_event[i][3][4] ) #dR(q2,q4)
    dR_patron_patron[i][5] = delta_R(quark_in_each_event[i][2][3], quark_in_each_event[i][2][4], quark_in_each_event[i][3][3], quark_in_each_event[i][3][4] #dR(q3,q4)




for i in range(0,10):
    dR_patron_jet[i][0] = delta_R(quark_in_each_event[i][0][3], quark_in_each_event[i][0][4], jet.eta[i][0], jet.phi[i][0] )
    dR_patron_jet[i][1] = delta_R(quark_in_each_event[i][1][3], quark_in_each_event[i][1][4], jet.eta[i][0], jet.phi[i][0] )
    dR_patron_jet[i][2] = delta_R(quark_in_each_event[i][2][3], quark_in_each_event[i][2][4], jet.eta[i][0], jet.phi[i][0] )
    dR_patron_jet[i][3] = delta_R(quark_in_each_event[i][3][3], quark_in_each_event[i][3][4], jet.eta[i][0], jet.phi[i][0] )


for i in range(0,10):
    for j in range(0,6):
        print(dR_patron_patron[i][j])




for i in range(0,10):
    for j in range(0,4):
        print(dR_patron_jet[i][j])




def min_delta_R(target_1, target_2):
    pass

