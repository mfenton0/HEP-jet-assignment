import uproot, sys, hdf5plugin
import pandas as pd 
import numpy as np 
from particle_properties_uproot import particle_properties  #import particle properties helper function from particle_properties.py
from jet_properties_uproot import jet_properties  #import jet properties helper function from jet_properties.py
import h5py


FILE_PATH = sys.argv[1]
STORE_PATH = sys.argv[2]

data  = uproot.open(FILE_PATH)['Delphes']
#data.show()

particle = particle_properties(data)
jet = jet_properties(data)

#Length = len(particle.event)
#test_length = 10

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

parton_array = np.zeros([ len(particle.event) , 6, 7])


#Generate maker for each stage(event selection and jet selection.)
marker_event = []
marker_jet = []

for i in range(len(particle.event)):
    marker_event.append(0)
    marker_jet.append(np.zeros([len(jet.pt[i])]))


marker_event = np.asanyarray(marker_event)
marker_jet = np.asanyarray(marker_jet)

print(type(marker_event), type(marker_jet))
print(marker_event.shape, marker_jet.shape)



#Mark which event pass the selection
print("+-----------------------------------------------------------------------------------------------------+")
print("Start event selection.")
for i in range(len(particle.event)):
    min_pt = np.min(jet.pt[i])
    num_of_eta_in_range = np.sum(jet.eta[i] < 2.4 ) 
    num_of_jet = len(jet.pt[i])
    num_of_btagged = np.sum(jet.btag[i] == 1)
    if min_pt > 20 and num_of_eta_in_range >= 6 and num_of_jet >=6 and num_of_btagged >= 2: 
        marker_event[i] = 1
    else :
        pass
print("Event selection doen.")
print("+-----------------------------------------------------------------------------------------------------+")

#Mark which jet in each event pass the selection.
print("+-----------------------------------------------------------------------------------------------------+")
print("Start jet selection.")
for i in range(len(particle.event)):
    if marker_event[i] == 1:
        for j in range(len(jet.pt[i])):
            if jet.btag[i][j] == 1 and jet.pt[i][j] > 20 and jet.eta[i][j] < 2.4:
                marker_jet[i][j] = 1 
            elif jet.pt[i][j] > 20 and jet.eta[i][j] <= 2.4:
                marker_jet[i][j] = 1
            else :
                pass
        else :
            pass 
print("Jet selection doen.")
print("+-----------------------------------------------------------------------------------------------------+")

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


for i in range(len(particle.event)):
    print("+------------------------------------------------------------------------------------------------------+")
    print("Start parsing event : {0}\nStart to trace top quark and find its daughters.".format(i))
    top_idx[i], top_daughter_idx_1[i], top_daughter_pid_1[i], top_daughter_idx_2[i], top_daughter_pid_2[i] = particle_tracing(particle.dataframelize(i), PID_TOP, 22)
    print("+------------------------------------------------------~-----------------------------------------------+")
    print("Start to find top_bar quark and its daughters.")
    top_bar_idx[i], top_bar_daughter_idx_1[i], top_bar_daughter_pid_1[i], top_bar_daughter_idx_2[i], top_bar_daughter_pid_2[i] = particle_tracing(particle.dataframelize(i), PID_TOP_BAR, 22)
    print("+------------------------------------------------------------------------------------------------------+")


### Tracing the daughter
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

    
    if daughter_1_pid == mother_pid or daughter_2_pid == mother_pid:

        init_idx = W_boson_idx
        daughter_pid = daughter_1_pid
        if daughter_2_pid == mother_pid:
            daughter_pid = daughter_2_pid
        while daughter_pid == mother_pid :
            daughter_1_idx = dataset.iloc[int(init_idx), 4]
            daughter_2_idx = dataset.iloc[int(init_idx), 5]

            daughter_1_pid = dataset.iloc[int(daughter_1_idx), 6]
            daughter_2_pid = dataset.iloc[int(daughter_2_idx), 6]

            daughter_pid = daughter_1_pid
            init_idx = daughter_1_idx
            if daughter_2_pid == mother_pid:
                daughter_pid = daughter_2_pid
                init_idx = daughter_2_idx
            
            
            print("Temporary daughter 1 indxe: {0}, PID: {1}".format(daughter_1_idx, daughter_1_pid))
            print("Temporary daughter 2 indxe: {0}, PID: {1}".format(daughter_2_idx, daughter_2_pid))

    
    print("Found daughter 1 index: {0}, PID: {1}.\nFound daughter 2 index: {2}, PID: {3}".format(daughter_1_idx, daughter_1_pid, daughter_2_idx, daughter_2_pid))
    return  b_quark_idx, daughter_1_idx, daughter_2_idx


for i in range(len(particle.event)):
    if marker_event[i] == 1 :
        print("+------------------------------------------------------------------------------------------------------+")
        print("Start parsing event : {0}\nStart to find top quark's daughters.".format(i))
        parton_array[i][0][0], parton_array[i][1][0], parton_array[i][2][0] = quark_finder(particle.dataframelize(i), top_daughter_idx_1[i], top_daughter_idx_2[i])
        print("+------------------------------------------------------~-----------------------------------------------+")
        print("Start to find top_bar quark's daughters.")
        parton_array[i][3][0], parton_array[i][4][0], parton_array[i][5][0], = quark_finder(particle.dataframelize(i), top_bar_daughter_idx_1[i], top_bar_daughter_idx_2[i])
        print("+------------------------------------------------------------------------------------------------------+")
    elif marker_event[i] == 0 :
        parton_array[i] = 'Nan'
    else: pass

 
# t t~ W+ W- b b~ 
# 0 0  0  0  0 0

# i.e.
# col 3 = 100010 
# col 5 = 101000
# col 6 = 101000
# col 8 = 010100
# col 10= 010001
# col 11= 010001

barcode = np.array([34, 40, 40, 17, 20, 20])
for i in range(len(particle.event)):
    if marker_event[i] == 1:
        for j in range(0,6):
            dataset = particle.dataframelize(i)
            parton_array[i][j][1] = dataset.iloc[int(parton_array[i][j][0]), 6]  #PDGID
            parton_array[i][j][2] = barcode[j]
            parton_array[i][j][3] = dataset.iloc[int(parton_array[i][j][0]), 7]  #Pt
            parton_array[i][j][4] = dataset.iloc[int(parton_array[i][j][0]), 8]  #Eta
            parton_array[i][j][5] = dataset.iloc[int(parton_array[i][j][0]), 9]  #Phi
            parton_array[i][j][6] = dataset.iloc[int(parton_array[i][j][0]), 10]  #Mass

def deltaPhi(phi1,phi2):
    phi = phi1-phi2
    while phi >= np.pi: phi -= np.pi*2.
    while phi < -np.pi: phi += np.pi*2.
    return phi

def delta_R(eta1, phi1, eta2, phi2):
    return np.sqrt(deltaPhi(phi1,phi2)**2+(eta1-eta2)**2)

def min_delta_R(target_1, target_2):
    pass

dR_between_parton_jet = []
dR_between_parton_parton = []

for i in range(len(particle.event)):
    dR_between_parton_jet.append(np.zeros([len(jet.pt[i]) * 6])) # # of connection = num of jet * num of parton
    dR_between_parton_parton.append(np.zeros([15])) # C^{6}_{2} = 15

dR_between_parton_jet = np.asanyarray(dR_between_parton_jet)
dR_between_parton_parton = np.asanyarray(dR_between_parton_parton)


max_num_of_jet_cand = []
for i in range(len(particle.event)):
    max_num_of_jet_cand.append(len(jet.pt[i]))
max_num_of_jet_cand = np.asanyarray(max_num_of_jet_cand)
max_num_of_jet = max_num_of_jet_cand.max()
print(max_num_of_jet)

#parton_jet_matching = np.zeros([len(jet.event), 6, 2])
matching_jet = []
matching_parton = []
for i in range(len(particle.event)):
    matching_jet.append(np.zeros([len(jet.pt[i])]))
    matching_parton.append(np.zeros([6]))

matching_jet = np.array(matching_jet)
matching_parton = np.array(matching_parton)


for i in range(len(particle.event)):
    if marker_event[i] == 1:
        j = 0
        a = 0
        b = 0
        while a < 6 :
            for b in range( len(jet.pt[i]) ):
                print(i, a, b)
                print(delta_R( parton_array[i][a][4], parton_array[i][a][5], jet.eta[i][b], jet.phi[i][b]))
                dR_between_parton_jet[i][j] = delta_R( parton_array[i][a][4], parton_array[i][a][5], jet.eta[i][b], jet.phi[i][b])
                j +=1
            a += 1 
    else :
        dR_between_parton_jet[i] = 'Nan'

for i in range(len(particle.event)):
    if marker_event[i] == 1:
        print("+------------------------------------------------------------------------------------------------------+")
        # print(dR_between_parton_jet.shape)
        array = np.reshape(dR_between_parton_jet[i], [6, len(jet.pt[i])])
        print(array.shape)
        
        dataset = pd.DataFrame({'0': array[0,:], 
                                '1': array[1,:],
                                '2': array[2,:],
                                '3': array[3,:],
                                '4': array[4,:],
                                '5': array[5,:],
                                })
        print(dataset)

        for j in range(0,6):
            print("+------------------------------------------------------------------------------------------------------+")
            min_val = dataset.stack().min()
            if min_val < 0.4:
                print("Min val: {0}".format(min_val))
                min_idx, min_col = dataset.stack().idxmin()
                matching_parton[i][j] = int(min_idx)
                matching_jet[i][j] = int(min_col)
                #parton_jet_matching[i][j][0] = int(min_idx)
                #parton_jet_matching[i][j][1] = int(min_col)
                print("The position of minimun appears. Raw: {0}, Colume: {1}".format(min_idx, min_col))
                dataset = dataset.drop([min_col], axis=1)
                dataset = dataset.drop([min_idx], axis=0)
                print("The dataset after delete the minimun's raw and colume:")
                print(dataset)
            else:
                matching_parton[i][j] = 'Nan'
                matching_jet[i][j] = 'Nan'
                #parton_jet_matching[i][j][0] = 'Nan'
                #parton_jet_matching[i][j][1] = 'Nan'
        for k in range(6, len(jet.pt[i])):
            matching_jet[i][k] = 'Nan'
    else : pass

parton_index = np.zeros([len(jet.event), 6])
jet_index = []
np.zeros([len(jet.event), 6])
for i in range(len(particle.event)):
    jet_index.append(np.zeros([len(jet.pt[i])]))


for i in range(len(particle.event)):
    if marker_event[i] == 1:
        for j in range(0,6):
            parton_index[i][j] = matching_parton[i][j]
        for k in range(len(jet.pt[i])):
            jet_index[i][k] = matching_jet[i][k]

jet_barcode = []
for i in range(len(particle.event)):
    jet_barcode.append(np.zeros([len(jet.pt[i])]))

jet_barcode = np.array(jet_barcode)

for i in range(len(particle.event)):
    if marker_event[i] == 1:
        for j in range(len(jet_index[i])):
            if jet_index[i][j] == 0:
                jet_barcode[i][j] = barcode[0]
            elif jet_index[i][j] == 1: 
                jet_barcode[i][j] = barcode[1]
            elif jet_index[i][j] == 2: 
                jet_barcode[i][j] = barcode[2]
            elif jet_index[i][j] == 3: 
                jet_barcode[i][j] = barcode[3]
            elif jet_index[i][j] == 4: 
                jet_barcode[i][j] = barcode[4]
            elif jet_index[i][j] == 5: 
                jet_barcode[i][j] = barcode[5]
            else :
                jet_barcode[i][j] = 'Nan'

jet_pt = []
jet_eta = []
jet_phi = []
jet_btag = []
jet_mass = []

for i in range(len(particle.event)):
    jet_pt.append(np.zeros([len(jet.pt[i])]))
    jet_eta.append(np.zeros([len(jet.pt[i])]))
    jet_phi.append(np.zeros([len(jet.pt[i])]))
    jet_btag.append(np.zeros([len(jet.pt[i])]))
    jet_mass.append(np.zeros([len(jet.pt[i])]))

jet_pt = np.array(jet_pt)
jet_eta = np.array(jet_eta)
jet_phi = np.array(jet_phi)
jet_btag = np.array(jet_btag)
jet_mass = np.array(jet_mass)

for i in range(len(particle.event)):
    if marker_event[i] == 1:
        for j in range(len(jet.pt[i])):
            if marker_jet[i][j] == 1:
                jet_pt[i][j] = jet.pt[i][j]
                jet_eta[i][j] = jet.eta[i][j]
                jet_phi[i][j] = jet.phi[i][j]
                jet_btag[i][j] = jet.btag[i][j]
                jet_mass[i][j] = jet.mass[i][j]
            else :
                jet_pt[i][j] = 'Nan'
                jet_eta[i][j] = 'Nan'
                jet_phi[i][j] = 'Nan'
                jet_btag[i][j] = 'Nan'
                jet_mass[i][j] = 'Nan'

hdf5_jet_parton_index = []
hdf5_jet_barcode = []
hdf5_jet_pt = []
hdf5_jet_eta = []
hdf5_jet_phi = []
hdf5_jet_mass = []
hdf5_jet_btagged = []

hdf5_parton_jet_index = []
hdf5_parton_pdgid = []
hdf5_parton_barcode = []
hdf5_parton_pt = []
hdf5_parton_eta = []
hdf5_parton_phi = []
hdf5_parton_mass = []


for i in range(len(particle.event)):
    if marker_event[i] == 1:
        hdf5_jet_parton_index.append(parton_index[i])
        hdf5_jet_barcode.append(jet_barcode[i])
        hdf5_jet_pt.append(jet_pt[i])
        hdf5_jet_eta.append(jet_eta[i])
        hdf5_jet_phi.append(jet_phi[i])
        hdf5_jet_mass.append(jet_mass[i])
        hdf5_jet_btagged.append(jet_btag[i])
    else: pass



for i in range(len(particle.event)):
    if marker_event[i] == 1:
        parton_pdgid = []
        parton_pt = []
        parton_eta = []
        parton_phi = []
        parton_mass = []
        for j in range(0,6):
            parton_pdgid.append(parton_array[i][j][1])
            parton_pt.append(parton_array[i][j][3])
            parton_eta.append(parton_array[i][j][4])
            parton_phi.append(parton_array[i][j][5])
            parton_mass.append(parton_array[i][j][6])

        hdf5_parton_jet_index.append(jet_index[i])
        hdf5_parton_pdgid.append(parton_pdgid)
        hdf5_parton_barcode.append(barcode)
        hdf5_parton_pt.append(parton_pt)
        hdf5_parton_eta.append(parton_eta)
        hdf5_parton_phi.append(parton_phi)
        hdf5_parton_mass.append(parton_mass)


hdf5_jet_parton_index = np.array(hdf5_jet_parton_index)
hdf5_jet_barcode = np.array(hdf5_jet_barcode)
hdf5_jet_pt = np.array(hdf5_jet_pt)
hdf5_jet_eta = np.array(hdf5_jet_eta)
hdf5_jet_phi = np.array(hdf5_jet_phi)
hdf5_jet_mass = np.array(hdf5_jet_mass)
hdf5_jet_btagged = np.array(hdf5_jet_btagged)

hdf5_parton_jet_index = np.array(hdf5_parton_jet_index)
hdf5_parton_pdgid = np.array(hdf5_parton_pdgid)
hdf5_parton_barcode = np.array(hdf5_parton_barcode)
hdf5_parton_pt = np.array(hdf5_parton_pt)
hdf5_parton_eta = np.array(hdf5_parton_eta)
hdf5_parton_phi = np.array(hdf5_parton_phi)
hdf5_parton_mass = np.array(hdf5_parton_mass)

lene = len(hdf5_jet_parton_index)


#Save the event which pass the selection
with h5py.File(STORE_PATH,'w') as f:

    dt = h5py.vlen_dtype(np.dtype('int32'))

    jet_parton_index = f.create_dataset('jet_parton_index', (lene, ), dtype=dt)
    jet_barcode = f.create_dataset('jet_barcode', (lene, ), dtype=dt)
    jet_pt = f.create_dataset('jet_pt', (lene, ), dtype=dt)
    jet_eta = f.create_dataset('jet_eta', (lene, ), dtype=dt)
    jet_phi = f.create_dataset('jet_phi', (lene, ), dtype=dt)
    jet_mass = f.create_dataset('jet_mass', (lene, ), dtype=dt)
    jet_btag = f.create_dataset('jet_btag', (lene, ), dtype=dt)

    for i in range(lene):
        jet_parton_index[i] = hdf5_jet_parton_index[i]
        jet_barcode[i] = hdf5_jet_parton_index[i]
        jet_pt[i] = hdf5_jet_pt[i]
        jet_eta[i] = hdf5_jet_eta[i]
        jet_phi[i] = hdf5_jet_phi[i]
        jet_mass[i] = hdf5_jet_mass[i]
        jet_btag[i] = hdf5_jet_btagged[i]

    parton_jet_index = f.create_dataset('hdf5_parton_jet_index', (lene, ), dtype=dt)
    parton_pdgid = f.create_dataset('hdf5_parton_pdgid', (lene, ), dtype=dt)
    parton_barcode = f.create_dataset('hdf5_parton_barcode', (lene, ), dtype=dt)
    parton_pt = f.create_dataset('hdf5_parton_pt', (lene, ), dtype=dt)
    parton_eta = f.create_dataset('hdf5_parton_eta', (lene, ), dtype=dt)
    parton_phi = f.create_dataset('hdf5_parton_phi', (lene, ), dtype=dt)
    parton_mass = f.create_dataset('hdf5_parton_mass', (lene, ), dtype=dt)

    for i in range(lene):
        parton_jet_index[i] = hdf5_parton_jet_index[i]
        parton_pdgid[i] = hdf5_parton_pdgid[i]
        parton_barcode[i] = hdf5_parton_barcode[i]
        parton_pt[i] = np.array(hdf5_parton_pt[i])
        parton_eta[i] = np.array(hdf5_parton_eta[i])
        parton_phi[i] = np.array(hdf5_parton_phi[i])
        parton_mass[i] = np.array(hdf5_parton_mass[i])