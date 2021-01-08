"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
#Import packages
import uproot
import pandas as pd 
import numpy as np 
from .particle_properties_uproot import particle_properties  #import particle properties helper function from particle_properties.py
from .jet_properties_uproot import jet_properties  #import jet properties helper function from jet_properties.py
import h5py, sys, traceback, os, tqdm
from .utilize import delta_R, deltaPhi, pdgid, event_selection, quark_finder, particle_tracing, deltaR_matching, barcode_recorder

class properties_for_particle():
    def __init__(self, event, pt, eta, phi, pid, M1, M2, D1, D2, status, rapidity, mass, charge):
        self.event = event
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.pid = pid
        self.M1 = M1
        self.M2 = M2
        self.D1 = D1
        self.D2 = D2
        self.status = status
        self.rapidity = rapidity
        self.mass = mass
        self.charge = charge
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

class properties_for_jet():
    def __init__(self, event, pt, eta, phi, btag, area, mass, charge):
        self.event = event
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.btag = btag
        self.area = area
        self.mass = mass
        self.charge = charge
    def dataframelize(self, index):

        idx = np.linspace(0, len( self.pt[index])-1, num = len( self.pt[index]) )

        jet_dict = {
                "Index": idx,
                "PT":  self.pt[index],
                "Eta":  self.eta[index],
                "Phi":  self.phi[index],
                "Mass":  self.mass[index],
                "Btag": self.btag[index],
                "Area": self.area[index]
            }
        jet_df = pd.DataFrame(jet_dict)
        return jet_df

def parse(INPUT_FILE, OUTPUT_FILE, MODEL, SINGLE):

    PID = pdgid()

    if int(SINGLE) == 1:
        try:
            print("Loading root file from {0}.".format(INPUT_FILE))
            data  = uproot.open(INPUT_FILE)['Delphes']
        except:
            print('Please check input file path.')
        particle = particle_properties(data)
        jet = jet_properties(data)

    elif int(SINGLE) != 1 and bool(SINGLE.isdigit) == True:
        _data = []
        files = os.listdir(INPUT_FILE)
        num_of_files = len(files)
        pbar = tqdm.tqdm(total=num_of_files)
        for a in files:
            try:
                print("Loading root file from {0}.".format(os.path.join(INPUT_FILE, a)))
                tmp = uproot.open(os.path.join(INPUT_FILE, a))['Delphes']
                _data.append(tmp)
                pbar.update(1) 
            except:
                print('Please check input file path.')
        
        for i in tqdm.trange(num_of_files):
            tmp_particle = particle_properties(_data[i])
            tmp_jet = jet_properties(_data[i])

            if i == 0 :
                _particle_event = tmp_particle.event
                _particle_pt = tmp_particle.pt
                _particle_eta = tmp_particle.eta
                _particle_phi = tmp_particle.phi
                _particle_pid = tmp_particle.pid
                _particle_M1 = tmp_particle.M1
                _particle_M2 = tmp_particle.M2
                _particle_D1 = tmp_particle.D1
                _particle_D2 = tmp_particle.D2
                _particle_status = tmp_particle.status
                _particle_rapidity = tmp_particle.rapidity
                _particle_mass = tmp_particle.mass
                _particle_charge = tmp_particle.charge

                _jet_event = tmp_jet.event
                _jet_pt = tmp_jet.pt
                _jet_eta = tmp_jet.eta
                _jet_phi = tmp_jet.phi
                _jet_btag = tmp_jet.btag
                _jet_area = tmp_jet.area
                _jet_mass = tmp_jet.mass
                _jet_charge = tmp_jet.charge
                
            else :
                _particle_event = np.concatenate((_particle_event, tmp_particle.event))
                _particle_pt = np.concatenate((_particle_pt, tmp_particle.pt))
                _particle_eta = np.concatenate((_particle_eta, tmp_particle.eta))
                _particle_phi = np.concatenate((_particle_phi, tmp_particle.phi))
                _particle_pid = np.concatenate((_particle_pid, tmp_particle.pid))
                _particle_M1 = np.concatenate((_particle_M1, tmp_particle.M1))
                _particle_M2 = np.concatenate((_particle_M2, tmp_particle.M2))
                _particle_D1 = np.concatenate((_particle_D1, tmp_particle.D1))
                _particle_D2 = np.concatenate((_particle_D2, tmp_particle.D2))
                _particle_status = np.concatenate((_particle_status, tmp_particle.status))
                _particle_rapidity = np.concatenate((_particle_rapidity, tmp_particle.rapidity))
                _particle_mass = np.concatenate((_particle_mass, tmp_particle.mass))
                _particle_charge = np.concatenate((_particle_charge, tmp_particle.charge))

                _jet_event = np.concatenate((_jet_event, tmp_jet.event))
                _jet_pt = np.concatenate((_jet_pt, tmp_jet.pt))
                _jet_eta = np.concatenate((_jet_eta, tmp_jet.eta))
                _jet_phi = np.concatenate((_jet_phi, tmp_jet.phi))
                _jet_btag = np.concatenate((_jet_btag, tmp_jet.btag))
                _jet_area = np.concatenate((_jet_area, tmp_jet.area))
                _jet_mass = np.concatenate((_jet_mass, tmp_jet.mass))
                _jet_charge = np.concatenate((_jet_charge, tmp_jet.charge))

        particle = properties_for_particle(_particle_event, _particle_pt, _particle_eta, _particle_phi, _particle_pid, _particle_M1, _particle_M2, _particle_D1, _particle_D2, _particle_status, _particle_rapidity, _particle_mass, _particle_charge)

        jet = properties_for_jet(_jet_event, _jet_pt, _jet_eta, _jet_phi, _jet_btag, _jet_area, _jet_mass, _jet_charge)
    else:
        print("Please inpurt a correct mode.\n1.-s 1\n2.-s > 1")

    
    if MODEL == "ttbar":
        """
        Barcode system
        t t~ W+ W- b b~ 
        0 0  0  0  0 0
        daughter of top and W+: 101000 ----> 40
        daughter of top and b: 101000 ----> 34
        daughter of anti top and W-: 100100 ----> 20
        daughter of anti top and b~: 100001 ----> 17
        """
        barcode = np.array([34, 40, 40, 17, 20, 20])
        NUM_OF_PARTON = 6
        
    elif MODEL == "ttH":
        """
        Barcode system
        t t~ W+ W- b b~ H
        0 0  0  0  0 0  0
        daughter of t and b = 1000100  ----> 68
        daughter of t and W+ = 1010000 ----> 80
        daughter of t~ and W- = 0101000 ----> 34
        daughter of t~ and b~ = 0100010 ----> 40
        daughter of H = 0000001 ----> 1
        """
        barcode = np.array([68, 80, 80, 34, 40, 40, 1, 1])
        NUM_OF_PARTON = 8
    elif MODEL == "four_top":
        """
        Barcode system
        t1 t2 t1~ t2~ W+1 W-1 W+2 W-2 b1 b2 b1~ b2~             describe          barcode   sequence
        0  0   0   0   0   0   0   0  0  0   0   0

        1  0   0   0   1   0   0   0  0  0   0   0  <--- daughter of t1 and W+1   2176         2,3
        1  0   0   0   0   0   0   0  1  0   0   0  <--- daughter of t1 and b1    2056          1
        0  0   1   0   0   1   0   0  0  0   0   0  <--- daughter of t1~ and W-1  576          5,6
        0  0   1   0   0   0   0   0  0  1   0   0  <--- daughter of t1~ and b1~  516           4

        0  1   0   0   0   0   1   0  0  0   0   0  <--- daughter of t2 and W+2   1056         7,8
        0  1   0   0   0   0   0   0  0  1   0   0  <--- daughter of t2 and b2    1028          9
        0  0   0   1   0   0   0   1  0  0   0   0  <--- daughter of t2~ and W-2  272          11.12
        0  0   0   1   0   0   0   0  0  0   0   1  <--- daughter of t2~ and b2~  257           10
        
        """
        barcode = np.array([2056, 2176, 2176, 516, 576, 576, 1028, 1056, 1056, 257, 272, 272])
        NUM_OF_PARTON = 12
        
    else:
        print("Please select a correct model.")

    #Mark which jet in each event pass the selection.
    print("+------------------------------------------------------------------------------------------------------+")
    print("Start jet selection.")
    print("+------------------------------------------------------------------------------------------------------+")
    marker_event, marker_jet, marker_btag = event_selection(jet.pt, jet.eta, jet.btag, MODEL)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Jet selection doen. {0} events has been selected.".format(np.sum(marker_event == 1)))
    print("+------------------------------------------------------------------------------------------------------+")

    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording the kinematics variables of jets in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    #Record the kinematical variables of jet in the selected event.
    jet_pt = []
    jet_eta = []
    jet_phi = []
    jet_mass = []
    jet_btag = []

    for i in tqdm.trange(len(jet.event)):
        if marker_event[i] == 1:
            jet_pt_tmp = []
            jet_eta_tmp = []
            jet_phi_tmp = []
            jet_mass_tmp = []
            jet_btag_tmp = []
            for j in range(len(jet.pt[i])):
                jet_pt_tmp.append(jet.pt[i][j])
                jet_eta_tmp.append(jet.eta[i][j])
                jet_phi_tmp.append(jet.phi[i][j])
                jet_mass_tmp.append(jet.mass[i][j])
                jet_btag_tmp.append(jet.btag[i][j])
                    
            jet_pt.append(jet_pt_tmp)
            jet_eta.append(jet_eta_tmp)
            jet_phi.append(jet_phi_tmp)
            jet_mass.append(jet_mass_tmp)
            jet_btag.append(jet_btag_tmp)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Finished to record the kinematics variables of jets in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    
    print("+------------------------------------------------------------------------------------------------------+")
    print("Starting parton tracing and looking for its daughter.")
    print("+------------------------------------------------------------------------------------------------------+")
    #Particle tracing and daughter finding section
    if MODEL == 'ttbar':
        top_idx, top_daughter_idx_1, top_daughter_pid_1, top_daughter_idx_2, top_daughter_pid_2 = [], [], [], [], []
        top_bar_idx, top_bar_daughter_idx_1, top_bar_daughter_pid_1, top_bar_daughter_idx_2, top_bar_daughter_pid_2 = [], [], [], [], []
        for i in tqdm.trange(len(particle.event)):
            _top_idx, _top_daughter_idx_1, _top_daughter_pid_1, _top_daughter_idx_2, _top_daughter_pid_2  = 0, 0, 0, 0, 0
            _top_bar_idx, _top_bar_daughter_idx_1, _top_bar_daughter_pid_1, _top_bar_daughter_idx_2, _top_bar_daughter_pid_2  = 0, 0, 0, 0, 0
            if marker_event[i] == 1:
                a = particle.dataframelize(i)
                _top_idx, _top_daughter_idx_1, _top_daughter_pid_1, _top_daughter_idx_2, _top_daughter_pid_2  = particle_tracing(particle.dataframelize(i), PID.top, 22, MODEL)
                _top_bar_idx, _top_bar_daughter_idx_1, _top_bar_daughter_pid_1, _top_bar_daughter_idx_2, _top_bar_daughter_pid_2  = particle_tracing(particle.dataframelize(i), PID.anti_top, 22, MODEL)
            else: 
                pass
            top_idx.append(_top_idx)
            top_daughter_idx_1.append(_top_daughter_idx_1)
            top_daughter_pid_1.append(_top_daughter_pid_1)
            top_daughter_idx_2.append(_top_daughter_idx_2)
            top_daughter_pid_2.append(_top_daughter_pid_2)
            top_bar_idx.append(_top_bar_idx)
            top_bar_daughter_idx_1.append(_top_bar_daughter_idx_1)
            top_bar_daughter_pid_1.append(_top_bar_daughter_pid_1)
            top_bar_daughter_idx_2.append(_top_bar_daughter_idx_2)
            top_bar_daughter_pid_2.append(_top_bar_daughter_pid_2)

        parton_array = np.zeros([ len(particle.event) , NUM_OF_PARTON, 7]) #7: index, pdgid, barcode, pt, eta, phi, mass

        for i in tqdm.trange(len(particle.event)):

            if marker_event[i] == 1 :
                parton_array[i][0][0], parton_array[i][1][0], parton_array[i][2][0] = quark_finder(particle.dataframelize(i), top_daughter_idx_1[i], top_daughter_idx_2[i])
                parton_array[i][3][0], parton_array[i][4][0], parton_array[i][5][0], = quark_finder(particle.dataframelize(i), top_bar_daughter_idx_1[i], top_bar_daughter_idx_2[i])
            elif marker_event[i] == 0 :
                parton_array[i] = 'Nan'
            else: pass
    elif MODEL == 'ttH':
        top_idx, top_daughter_idx_1, top_daughter_pid_1, top_daughter_idx_2, top_daughter_pid_2 = [], [], [], [], []
        top_bar_idx, top_bar_daughter_idx_1, top_bar_daughter_pid_1, top_bar_daughter_idx_2, top_bar_daughter_pid_2 = [], [], [], [], []
        higgs_idx, higgs_daughter_idx_1, higgs_daughter_pid_1, higgs_daughter_idx_2, higgs_daughter_pid_2 = [], [], [], [], []
        for i in tqdm.trange(len(particle.event)):
            _top_idx, _top_daughter_idx_1, _top_daughter_pid_1, _top_daughter_idx_2, _top_daughter_pid_2  = 0, 0, 0, 0, 0
            _top_bar_idx, _top_bar_daughter_idx_1, _top_bar_daughter_pid_1, _top_bar_daughter_idx_2, _top_bar_daughter_pid_2  = 0, 0, 0, 0, 0
            _higgs_idx, _higgs_daughter_idx_1, _higgs_daughter_pid_1, _higgs_daughter_idx_2, _higgs_daughter_pid_2 = 0, 0, 0, 0, 0
            if marker_event[i] == 1:
                _top_idx, _top_daughter_idx_1, _top_daughter_pid_1, _top_daughter_idx_2, _top_daughter_pid_2  = particle_tracing(particle.dataframelize(i), PID.top, 22, MODEL)
                _top_bar_idx, _top_bar_daughter_idx_1, _top_bar_daughter_pid_1, _top_bar_daughter_idx_2, _top_bar_daughter_pid_2  = particle_tracing(particle.dataframelize(i), PID.anti_top, 22, MODEL)
                _higgs_idx, _higgs_daughter_idx_1, _higgs_daughter_pid_1, _higgs_daughter_idx_2, _higgs_daughter_pid_2 = particle_tracing(particle.dataframelize(i), PID.higgs, 22, MODEL)
            else: 
                pass
            top_idx.append(_top_idx)
            top_daughter_idx_1.append(_top_daughter_idx_1)
            top_daughter_pid_1.append(_top_daughter_pid_1)
            top_daughter_idx_2.append(_top_daughter_idx_2)
            top_daughter_pid_2.append(_top_daughter_pid_2)
            top_bar_idx.append(_top_bar_idx)
            top_bar_daughter_idx_1.append(_top_bar_daughter_idx_1)
            top_bar_daughter_pid_1.append(_top_bar_daughter_pid_1)
            top_bar_daughter_idx_2.append(_top_bar_daughter_idx_2)
            top_bar_daughter_pid_2.append(_top_bar_daughter_pid_2)
            higgs_idx.append(_higgs_idx)
            higgs_daughter_idx_1.append(_higgs_daughter_idx_1)
            higgs_daughter_pid_1.append(_higgs_daughter_pid_1)
            higgs_daughter_idx_2.append(_higgs_daughter_idx_2)
            higgs_daughter_pid_2.append(_higgs_daughter_pid_2)

        parton_array = np.zeros([ len(particle.event) , NUM_OF_PARTON, 7]) #7: index, pdgid, barcode, pt, eta, phi, mass

        for i in tqdm.trange(len(particle.event)):

            if marker_event[i] == 1 :
                parton_array[i][0][0], parton_array[i][1][0], parton_array[i][2][0] = quark_finder(particle.dataframelize(i), top_daughter_idx_1[i], top_daughter_idx_2[i])
                parton_array[i][3][0], parton_array[i][4][0], parton_array[i][5][0], = quark_finder(particle.dataframelize(i), top_bar_daughter_idx_1[i], top_bar_daughter_idx_2[i])
                parton_array[i][6][0], parton_array[i][7][0] = higgs_daughter_idx_1[i], higgs_daughter_idx_2[i]
            elif marker_event[i] == 0 :
                parton_array[i] = 'Nan'
            else: pass
    elif MODEL == 'four_top':
        top_1_idx, top_1_daughter_idx_1, top_1_daughter_pid_1, top_1_daughter_idx_2, top_1_daughter_pid_2 = [], [], [], [], []
        top_2_idx, top_2_daughter_idx_1, top_2_daughter_pid_1, top_2_daughter_idx_2, top_2_daughter_pid_2 = [], [], [], [], []
        top_1_bar_idx, top_1_bar_daughter_idx_1, top_1_bar_daughter_pid_1, top_1_bar_daughter_idx_2, top_1_bar_daughter_pid_2 = [], [], [], [], []
        top_2_bar_idx, top_2_bar_daughter_idx_1, top_2_bar_daughter_pid_1, top_2_bar_daughter_idx_2, top_2_bar_daughter_pid_2 = [], [], [], [], []

        for i in tqdm.trange(len(particle.event)):

            _top_1_idx, _top_2_idx, _top_1_daughter_idx_1, _top_1_daughter_pid_1, _top_1_daughter_idx_2, _top_1_daughter_pid_2,  _top_2_daughter_idx_1, _top_2_daughter_pid_1, _top_2_daughter_idx_2, _top_2_daughter_pid_2 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            _top_1_bar_idx, _top_2_bar_idx, _top_1_bar_daughter_idx_1, _top_1_bar_daughter_pid_1, _top_1_bar_daughter_idx_2, _top_1_bar_daughter_pid_2,  _top_2_bar_daughter_idx_1, _top_2_bar_daughter_pid_1, _top_2_bar_daughter_idx_2, _top_2_bar_daughter_pid_2  = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

            if marker_event[i] == 1:

                _top_1_idx, _top_2_idx, _top_1_daughter_idx_1, _top_1_daughter_pid_1, _top_1_daughter_idx_2, _top_1_daughter_pid_2,  _top_2_daughter_idx_1, _top_2_daughter_pid_1, _top_2_daughter_idx_2, _top_2_daughter_pid_2  = particle_tracing(particle.dataframelize(i), PID.top, 22, MODEL)

                _top_1_bar_idx, _top_2_bar_idx, _top_1_bar_daughter_idx_1, _top_1_bar_daughter_pid_1, _top_1_bar_daughter_idx_2, _top_1_bar_daughter_pid_2,  _top_2_bar_daughter_idx_1, _top_2_bar_daughter_pid_1, _top_2_bar_daughter_idx_2, _top_2_bar_daughter_pid_2 = particle_tracing(particle.dataframelize(i), PID.anti_top, 22, MODEL)
            else: 
                pass
            top_1_idx.append(_top_1_idx)
            top_1_daughter_idx_1.append(_top_1_daughter_idx_1)
            top_1_daughter_pid_1.append(_top_1_daughter_pid_1) 
            top_1_daughter_idx_2.append(_top_1_daughter_idx_2)
            top_1_daughter_pid_2.append(_top_1_daughter_pid_2)
            top_2_idx.append(_top_2_idx)
            top_2_daughter_idx_1.append(_top_2_daughter_idx_1) 
            top_2_daughter_pid_1.append(_top_2_daughter_pid_1)
            top_2_daughter_idx_2.append(_top_2_daughter_idx_2)
            top_2_daughter_pid_2.append(_top_2_daughter_pid_2)
            top_1_bar_idx.append(_top_1_bar_idx)
            top_1_bar_daughter_idx_1.append(_top_1_bar_daughter_idx_1) 
            top_1_bar_daughter_pid_1.append(_top_1_bar_daughter_pid_1)
            top_1_bar_daughter_idx_2.append(_top_1_bar_daughter_idx_2)
            top_1_bar_daughter_pid_2.append(_top_1_bar_daughter_pid_2)
            top_2_bar_idx.append(_top_2_bar_idx)
            top_2_bar_daughter_idx_1.append(_top_2_bar_daughter_idx_1) 
            top_2_bar_daughter_pid_1.append(_top_2_bar_daughter_pid_1)
            top_2_bar_daughter_idx_2.append(_top_2_bar_daughter_idx_2)
            top_2_bar_daughter_pid_2.append(_top_2_bar_daughter_pid_2)


        parton_array = np.zeros([ len(particle.event) , NUM_OF_PARTON, 7]) #7: index, pdgid, barcode, pt, eta, phi, mass

        for i in tqdm.trange(len(particle.event)):

            if marker_event[i] == 1 :
                parton_array[i][0][0], parton_array[i][1][0], parton_array[i][2][0] = quark_finder(particle.dataframelize(i), top_1_daughter_idx_1[i], top_1_daughter_idx_2[i])
                parton_array[i][3][0], parton_array[i][4][0], parton_array[i][5][0], = quark_finder(particle.dataframelize(i), top_2_daughter_idx_1[i], top_2_daughter_idx_2[i])
                parton_array[i][6][0], parton_array[i][7][0], parton_array[i][8][0], = quark_finder(particle.dataframelize(i), top_1_bar_daughter_idx_1[i], top_1_bar_daughter_idx_2[i])
                parton_array[i][9][0], parton_array[i][10][0], parton_array[i][11][0], = quark_finder(particle.dataframelize(i), top_2_bar_daughter_idx_1[i], top_2_bar_daughter_idx_2[i])
            elif marker_event[i] == 0 :
                parton_array[i] = 'Nan'
            else: pass
    else :
        print("Please select a correct model.")


    print("+------------------------------------------------------------------------------------------------------+")
    print("Parton tracing section complete. The daughter of W+/W- and bbbar has been found.")
    print("+------------------------------------------------------------------------------------------------------+")

    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording the kinematics variables of partons in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    parton_pdgid = []
    parton_barcode = []
    parton_pt = []
    parton_eta = []
    parton_phi = []
    parton_mass = []


    for i in tqdm.trange(len(particle.event)):
        if marker_event[i] == 1:
            _parton_pdgid = []
            _parton_barcode = []
            _parton_pt = []
            _parton_eta = []
            _parton_phi = []
            _parton_mass = []
            for j in range(NUM_OF_PARTON):
                dataset = particle.dataframelize(i)

                _parton_pdgid.append(dataset.iloc[int(parton_array[i][j][0]), 6])
                _parton_barcode.append(barcode[j])
                _parton_pt.append(dataset.iloc[int(parton_array[i][j][0]), 7])
                _parton_eta.append(dataset.iloc[int(parton_array[i][j][0]), 8])
                _parton_phi.append(dataset.iloc[int(parton_array[i][j][0]), 9])
                _parton_mass.append(dataset.iloc[int(parton_array[i][j][0]), 10])


            parton_pdgid.append(_parton_pdgid)
            parton_barcode.append(_parton_barcode)
            parton_pt.append(_parton_pt)
            parton_eta.append(_parton_eta)
            parton_phi.append(_parton_phi)
            parton_mass.append(_parton_mass)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Finished to record the kinematics variables of partons in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")

    print("+------------------------------------------------------------------------------------------------------+")
    print("Starting parton-jet matching.")
    print("+------------------------------------------------------------------------------------------------------+")
    jet_parton_index = []
    parton_jet_index = []
    for i in tqdm.trange(len(jet_pt)):
        _jet_parton_index, _parton_jet_index = deltaR_matching(NUM_OF_PARTON, len(jet_pt[i]), parton_eta[i], parton_phi[i], jet_eta[i], jet_phi[i], 0.4, MODEL)
        jet_parton_index.append(_jet_parton_index)
        parton_jet_index.append(_parton_jet_index)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Parton-jet matching finished.")
    print("+------------------------------------------------------------------------------------------------------+")
    jet_parton_index = np.asanyarray(jet_parton_index, dtype=object)
    parton_jet_index = np.asanyarray(parton_jet_index, dtype=object)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording barcode information.")
    print("+------------------------------------------------------------------------------------------------------+")

    jet_barcode = []
    for i in tqdm.trange(len(jet_pt)):
        _jet_barcode = barcode_recorder(jet_parton_index[i], MODEL)
        jet_barcode.append(_jet_barcode)

    print("+------------------------------------------------------------------------------------------------------+")
    print("Barcode information has beed record.")
    print("+------------------------------------------------------------------------------------------------------+")

    if MODEL == 'ttH':
        target = [i for i in range(6)]
        N_match_top_in_event = np.zeros([len(jet_pt)])
        for i in tqdm.trange(len(jet_parton_index)):
        
            intersetion = set(target).intersection(jet_parton_index[i])
            N_match_top_in_event[i] = len(intersetion) // 3

        N_match_higgs_in_event = np.zeros([len(jet_pt)])
        for i in range(len(jet_parton_index)):
            if 7 in jet_parton_index[i]:
                if 6 in jet_parton_index[i]:
                    N_match_higgs_in_event[i] = 1
        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2)  ))
        print("+------------------------------------------------------------------------------------------------------+")
    elif MODEL == 'ttbar':
        target = [i for i in range(NUM_OF_PARTON)]
        N_match_top_in_event = np.zeros([len(jet_pt)])
        for i in tqdm.trange(len(jet_parton_index)):
        
            intersetion = set(target).intersection(jet_parton_index[i])
            N_match_top_in_event[i] = len(intersetion) // 3
        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2)  ))
        print("+------------------------------------------------------------------------------------------------------+")
    elif MODEL == 'four_top':

        target = [i for i in range(NUM_OF_PARTON)]
        N_match_top_in_event = np.zeros([len(jet_pt)])
        for i in tqdm.trange(len(jet_parton_index)):
            
            intersetion = set(target).intersection(jet_parton_index[i])
            N_match_top_in_event[i] = len(intersetion) // 3
        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.\nFound {2} events with 3 ttbar candidate exist.\nFound {3} events with 4 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2), np.sum(N_match_top_in_event == 3), np.sum(N_match_top_in_event == 4)  ))
        print("+------------------------------------------------------------------------------------------------------+")
    else : pass
    
   
    # Save selected events

    lene = len(parton_pdgid)

    print("+------------------------------------------------------------------------------------------------------+")
    print("Writing event record to the hdf5 file.")
    print("+------------------------------------------------------------------------------------------------------+")

    # Save the event which pass the selection
    with h5py.File(OUTPUT_FILE,'w') as f:
        dt = h5py.vlen_dtype(np.dtype('float16'))
        
        hdf5_jet_parton_index = f.create_dataset('jet_parton_index', (lene, ), dtype=dt)
        hdf5_jet_barcode = f.create_dataset('jet_barcode', (lene, ), dtype=dt)
        hdf5_jet_pt = f.create_dataset('jet_pt', (lene, ), dtype=dt)
        hdf5_jet_eta = f.create_dataset('jet_eta', (lene, ), dtype=dt)
        hdf5_jet_phi = f.create_dataset('jet_phi', (lene, ), dtype=dt)
        hdf5_jet_mass = f.create_dataset('jet_mass', (lene, ), dtype=dt)
        hdf5_jet_btag = f.create_dataset('jet_btag', (lene, ), dtype=dt)
        print("Writing jet information to {0}".format(OUTPUT_FILE))
        for i in tqdm.trange(lene):
            hdf5_jet_parton_index[i] = jet_parton_index[i]
            hdf5_jet_barcode[i] = jet_barcode[i]
            hdf5_jet_pt[i] = jet_pt[i]
            hdf5_jet_eta[i] = jet_eta[i]
            hdf5_jet_phi[i] = jet_phi[i]
            hdf5_jet_mass[i] = jet_mass[i]
            hdf5_jet_btag[i] = jet_btag[i]

        hdf5_parton_jet_index = f.create_dataset('parton_jet_index', (lene, ), dtype=dt)
        hdf5_parton_pdgid = f.create_dataset('parton_pdgid', (lene, ), dtype=dt)
        hdf5_parton_barcode = f.create_dataset('parton_barcode', (lene, ), dtype=dt)
        hdf5_parton_pt = f.create_dataset('parton_pt', (lene, ), dtype=dt)
        hdf5_parton_eta = f.create_dataset('parton_eta', (lene, ), dtype=dt)
        hdf5_parton_phi = f.create_dataset('parton_phi', (lene, ), dtype=dt)
        hdf5_parton_mass = f.create_dataset('parton_mass', (lene, ), dtype=dt)
        hdf5_N_match_top_in_event = f.create_dataset('N_match_top_in_event', data = N_match_top_in_event)
        if MODEL == "ttH":
            N_match_higgs_in_event = f.create_dataset('N_match_higgs_in_event', data = N_match_higgs_in_event)
        
        print("Writing parton information to {0}".format(OUTPUT_FILE))
        for i in tqdm.trange(lene):
            hdf5_parton_jet_index[i] = parton_jet_index[i]
            hdf5_parton_pdgid[i] = parton_pdgid[i]
            hdf5_parton_barcode[i] = parton_barcode[i]
            hdf5_parton_pt[i] = parton_pt[i]
            hdf5_parton_eta[i] = parton_eta[i]
            hdf5_parton_phi[i] = parton_phi[i]
            hdf5_parton_mass[i] = parton_mass[i]

    print("+------------------------------------------------------------------------------------------------------+")
    print("Event record has been send to {0}.".format(OUTPUT_FILE))
    print("+------------------------------------------------------------------------------------------------------+")
