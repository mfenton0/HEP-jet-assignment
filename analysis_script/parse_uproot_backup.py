#!/usr/bin/env python
# coding: utf-8

# In[53]:


import uproot
import pandas as pd 
import numpy as np 
from particle_properties import particle_properties
from jet_properties import jet_properties



data  = uproot.open('./tag_1_delphes_events.root')['Delphes']
#data.show()

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

W_plus_idx = np.zeros(len(particle.event))

W_minus_idx = np.zeros(len(particle.event))

b_quark = np.zeros(len(particle.event))

b_bar_quark = np.zeros(len(particle.event))

quark_idx_1 = np.zeros(len(particle.event))

quark_idx_2 = np.zeros(len(particle.event))

quark_idx_3 = np.zeros(len(particle.event))

quark_idx_4 = np.zeros(len(particle.event))


for i in range(0,10):
    top_idx[i], top_daughter_idx_1[i], top_daughter_pid_1[i], top_daughter_idx_2[i], top_daughter_pid_2[i] = particle_tracing(particle.dataframelize(i), PID_TOP, 22)
    top_bar_idx[i], top_bar_daughter_idx_1[i], top_bar_daughter_pid_1[i], top_bar_daughter_idx_2[i], top_bar_daughter_pid_2[i] = particle_tracing(particle.dataframelize(i), PID_TOP_BAR, 22)



#Input two daughter of top/top_bar and find their daughter
def quark_finder(dataset, mother_idx_1, mother_idx_2):
    
    #Specific two daughter of top
    def W_b_specifier(dataset, input_1_idx, input_2_idx):
        if dataset.iloc[int(input_1_idx),6] == PID_W_plus or dataset.iloc[int(input_1_idx),6] == PID_W_minus :
            return int(input_1_idx), int(dataset.iloc[int(input_1_idx),6]), int(input_2_idx)
        elif dataset.iloc[int(input_1_idx),6] == PID_BOTTOM or dataset.iloc[int(input_1_idx),6] == PID_BOTTOM_BAR :
            return  int(input_2_idx), int(dataset.iloc[int(input_1_idx),6]), int(input_1_idx)
        else :
            print("Please check your data.")
    
    W_boson_idx, mother_pid, b_quark_idx = W_b_specifier(dataset, mother_idx_1, mother_idx_2)
    
    #Find the two daughters of boson
    
    daughter_1_idx = dataset.iloc[W_boson_idx, 4]
    daughter_1_pid = dataset.iloc[daughter_1_idx, 6]
    daughter_2_idx = dataset.iloc[W_boson_idx, 5]
    daughter_2_pid = dataset.iloc[daughter_2_idx, 6]

    
    if daughter_1_pid == mother_pid and daughter_2_pid == mother_pid:
        init_idx = W_boson_idx
        while daughter_1_pid == mother_pid:
            daughter_1_idx = dataset.iloc[int(init_idx), 4]
            daughter_1_pid = dataset.iloc[int(daughter_1_idx), 6]
            init_idx = daughter_1_idx
            print("Temporary daughter 1 indxe: {0}, PID: {1}".format(daughter_1_idx, daughter_1_pid))
        init_idx = W_boson_idx
        while daughter_2_pid == mother_pid:
            daughter_2_idx = dataset.iloc[int(init_idx), 5]
            daughter_2_pid = dataset.iloc[int(daughter_2_idx), 6]
            init_idx = daughter_2_idx
            print("Temporary daughter 2 indxe: {0}, PID: {1}".format(daughter_2_idx, daughter_2_pid))
    
    print("Found daughter 1 index: {0}, PID: {1}.\nFound daughter 2 index: {2}, PID: {3}".format(daughter_1_idx, daughter_1_pid, daughter_2_idx, daughter_2_pid))
    return W_boson_idx, b_quark_idx, daughter_1_idx, daughter_2_idx




# In[54]:


for i in range(0,10):
    print("+-----------------------------------------------------------------------------------------------------+")
    print("Start parsing event : {0}\nStart to find top quark's daughters.")
    W_plus_idx[i], b_quark[i], quark_idx_1[i], quark_idx_2[i] = quark_finder(particle.dataframelize(i), top_daughter_idx_1[i], top_daughter_idx_2[i])
    print("+-----------------------------------------------------------------------------------------------------+")
    print("Start to find top_bar quark's daughters.")
    W_minus_idx[i], b_bar_quark[i], quark_idx_3[i], quark_idx_4[i] = quark_finder(particle.dataframelize(i), top_bar_daughter_idx_1[i], top_bar_daughter_idx_2[i])
    print("+-----------------------------------------------------------------------------------------------------+")


# In[ ]:


df9 = particle.dataframelize(9)
df9.iloc[841,:]


# In[6]:


df9.iloc[855,:]


# In[34]:


df9.iloc[854,:]


# In[19]:


df9.iloc[856,:]


# In[32]:


df9.iloc[857,:]


# In[ ]:




