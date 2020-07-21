import uproot
import pandas as pd 
import numpy as np 

class particle_properties():
    def __init__(self, data):
        self.event = data.array('Event')
        self.pt = data.array('Particle.PT')
        self.eta = data.array('Particle.Eta')
        self.phi = data.array('Particle.Phi')
        self.pid = data.array('Particle.PID')
        self.M1 = data.array('Particle.M1')
        self.M2 = data.array('Particle.M2')
        self.D1 = data.array('Particle.D1')
        self.D2 = data.array('Particle.D2')
        self.status = data.array('Particle.Status')
        self.rapidity = data.array('Particle.Rapidity')
        self.mass = data.array('Particle.Mass')
        self.charge = data.array('Particle.Charge')
        print(len(self.pt[0]), len(self.eta[0]), len(self.M1[0]))
    def dataframelize(self, index):

        idx = np.linspace(0, len( self.pt[index])-1, num = len( self.pt[index]) )

        patron_dict = {
                "Index": idx,
                "Status":  self.status[index],
                "Mother_1":  self.M1[index],
                "Mother_2":  self.M2[index],
                "Daughter_1":  self.D1[index],
                "Daughter_2":  self.D2[index],
                "PID":  self.pid[index],
                "PT":  self.pt[index],
                "Eta":  self.eta[index],
                "Phi":  self.phi[index],
                "Mass":  self.mass[index]
            }
        patron_df = pd.DataFrame(patron_dict)
        return patron_df

particle = particle_properties(data)

df0 = particle.dataframelize(0)
df1 = particle.dataframelize(1)
df2 = particle.dataframelize(2)


class jet_properties():
    def __init__(self, data):
        self.event = data.array('Event')
        self.pt = data.array('Jet.PT')
        self.eta = data.array('Jet.Eta')
        self.phi = data.array('Jet.Phi')
        self.btag = data.array('Jet.BTag')
        self.area = data.array('Jet.Area')
        self.mass = data.array('Jet.Mass')
        self.charge = data.array('Jet.Charge')
jet = jet_properties(data)

print(len(jet.event), len(jet.pt), len(jet.pt[0]))

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

for i in range(0,10):
    top_idx[i], top_daughter_idx_1[i], top_daughter_pid_1[i], top_daughter_idx_2[i], top_daughter_pid_2[i] = particle_tracing(particle.dataframelize(i), PID_TOP, 22)
    top_bar_idx[i], top_bar_daughter_idx_1[i], top_bar_daughter_pid_1[i], top_bar_daughter_idx_2[i], top_bar_daughter_pid_2[i] = particle_tracing(particle.dataframelize(i), PID_TOP_BAR, 22)


quark_idx_1 = np.zeros(len(particle.event))

quark_idx_2 = np.zeros(len(particle.event))

quark_idx_3 = np.zeros(len(particle.event))

quark_idx_4 = np.zeros(len(particle.event))

def quark_finder(dataset, mother_idx):

    print("Mother index: {0}".format(mother_idx))
    
    

    daughter_idx_1 = dataset.iloc[int(mother_idx), 4]
    daughter_idx_2 = dataset.iloc[int(mother_idx), 5]

    status = dataset.iloc[int(daughter_idx_1),1]
    print("Daughter's of t/t~'s status code: {0}".format(status))

    while status != 23 :
        mother_idx = daughter_idx_1
        daughter_idx_1 = dataset.iloc[int(mother_idx), 4]
        daughter_idx_2 = dataset.iloc[int(mother_idx), 5]
        status = dataset.iloc[int(daughter_idx_1), 1]
        #print("2nd stage daughter status code: {0}".format(status))

    
    daughter_pid_1 = dataset.iloc[daughter_idx_1, 6]
    daughter_pid_2 = dataset.iloc[daughter_idx_2, 6]
    
    print("Daughter 1 index: {0}, daughter 2 index: {1}".format(daughter_idx_1, daughter_idx_2))
    print("Daughter 1 PID: {0}, daughter 2 PID: {1}".format(daughter_pid_1, daughter_pid_2))
    
    return int(daughter_idx_1), int(daughter_idx_2)

for i in range(0,10):
    quark_idx_1[i], quark_idx_2[i] = quark_finder(particle.dataframelize(i), top_daughter_idx_1[i])
    quark_idx_3[i], quark_idx_4[i] = quark_finder(particle.dataframelize(i), top_bar_daughter_idx_1[i])




quark_in_each_event = np.zeros([len(particle.event), 4, 6])
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


