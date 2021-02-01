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
from .MissingET import Missing_ET
from .lepton import Lepton
import h5py, sys, traceback, os, tqdm, time
from .utilize import delta_R, deltaPhi, pdgid, event_selection, quark_finder, particle_tracing, deltaR_matching, barcode_recorder, deltaPhi
import multiprocessing as mp

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

class property_for_lepton():
    def __init__(self, e_pt, e_eta, e_phi, m_pt, m_eta, m_phi):
        self.electron_pt = e_pt
        self.electron_eta = e_eta
        self.electron_phi = e_phi
        self.Muon_pt = m_pt
        self.Muon_eta = m_eta
        self.Muon_phi = m_phi
        
    def dataframelize(self, index):

        idx = np.linspace(0, len( self.electron_pt[index])-1, num = len( self.electron_pt[index]) )

        Lepton_dict = {
                "Index": idx,
                "Electron_pt":  self.electron_pt[index],
                "Electron_eta":  self.electron_eta[index],
                "Electron_phi":  self.electron_phi[index],
                "Muon_pt":  self.Muon_pt[index],
                "Muon_eta":  self.Muon_eta[index],
                "Muon_phi":  self.Muon_phi[index],
            }
        Lepton_df = pd.DataFrame(Lepton_dict)
        return Lepton_df

class property_for_Missing_ET():
    def __init__(self, met, eta, phi):
        self.MET = met
        self.eta = eta
        self.phi = phi
        
    def dataframelize(self, index):

        idx = np.linspace(0, len( self.MET[index])-1, num = len( self.MET[index]) )

        MissingET_dict = {
                "Index": idx,
                "MET":  self.MET[index],
                "Eta":  self.eta[index],
                "Phi":  self.phi[index],

            }
        MissingET_df = pd.DataFrame(MissingET_dict)
        return MissingET_df

def parse(INPUT_FILE, OUTPUT_FILE, MODEL, SINGLE, PROCESS):

    PID = pdgid()

    if int(SINGLE) == 1:
        try:
            print("Loading root file from {0}.".format(INPUT_FILE))
            data  = uproot.open(INPUT_FILE)['Delphes']
        except:
            print('Please check input file path.')
        particle = particle_properties(data)
        jet = jet_properties(data)
        lepton = Lepton(data)
        missing_et = Missing_ET(data)

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
            tmp_lepton = Lepton(_data[i])
            tmp_missing_et = Missing_ET(_data[i])

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

                _lepton_electron_pt = tmp_lepton.electron_pt
                _lepton_electron_eta = tmp_lepton.electron_eta
                _lepton_electron_phi = tmp_lepton.electron_phi
                _lepton_muon_pt = tmp_lepton.Muon_pt
                _lepton_muon_eta = tmp_lepton.Muon_eta
                _lepton_muon_phi = tmp_lepton.Muon_phi

                _missing_et_met = tmp_missing_et.MET
                _missing_et_eta = tmp_missing_et.eta
                _missing_et_phi = tmp_missing_et.phi
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
                
                _lepton_electron_pt = np.concatenate((_lepton_electron_pt, tmp_lepton.electron_pt))
                _lepton_electron_eta = np.concatenate((_lepton_electron_eta, tmp_lepton.electron_eta))
                _lepton_electron_phi = np.concatenate((_lepton_electron_phi, tmp_lepton.electron_phi))
                _lepton_muon_pt = np.concatenate((_lepton_muon_pt, tmp_lepton.Muon_pt))
                _lepton_muon_eta = np.concatenate((_lepton_muon_eta, tmp_lepton.Muon_eta))
                _lepton_muon_phi = np.concatenate((_lepton_muon_phi, tmp_lepton.Muon_phi))

                _missing_et_met = np.concatenate((_missing_et_met, tmp_missing_et.MET))
                _missing_et_eta = np.concatenate((_missing_et_eta, tmp_missing_et.eta))
                _missing_et_phi = np.concatenate((_missing_et_phi, tmp_missing_et.phi))

        particle = properties_for_particle(_particle_event, _particle_pt, _particle_eta, _particle_phi, _particle_pid, _particle_M1, _particle_M2, _particle_D1, _particle_D2, _particle_status, _particle_rapidity, _particle_mass, _particle_charge)
        jet = properties_for_jet(_jet_event, _jet_pt, _jet_eta, _jet_phi, _jet_btag, _jet_area, _jet_mass, _jet_charge)
        lepton = property_for_lepton(_lepton_electron_pt, _lepton_electron_eta, _lepton_electron_phi, _lepton_muon_pt, _lepton_muon_eta, _lepton_muon_phi)
        missing_et = property_for_Missing_ET(_missing_et_met, _missing_et_eta, _missing_et_phi)
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
        NUM_OF_DAUGHTER = 6
    elif MODEL == "ttbar_lep_left":
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
        NUM_OF_PARTON = 4
        NUM_OF_DAUGHTER = 6
    elif MODEL == "ttbar_lep_right":
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
        NUM_OF_PARTON = 4
        NUM_OF_DAUGHTER = 6
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
        NUM_OF_DAUGHTER = 8
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
        NUM_OF_DAUGHTER = 12
    else:
        print("Please select a correct model.")

    tmp_lepton_pt = []
    tmp_lepton_eta = []
    tmp_lepton_phi = []
    tmp_lepton_pdgid = []

    for i in range(len(lepton.electron_pt)):
        _lepton_pt = []
        _lepton_eta = []
        _lepton_phi = []
        _lepton_pdgid = []
        if len(lepton.electron_pt[i]) != 0:
            for a, b, c in zip(lepton.electron_pt[i], lepton.electron_eta[i], lepton.electron_phi[i]):
                _lepton_pt.append(a)
                _lepton_eta.append(b)
                _lepton_phi.append(c)
                _lepton_pdgid.append(11)
        else:
            pass
        if len(lepton.Muon_pt[i]) != 0:
            for a, b, c in zip(lepton.Muon_pt[i], lepton.Muon_eta[i], lepton.Muon_phi[i]):
                _lepton_pt.append(a)
                _lepton_eta.append(b)
                _lepton_phi.append(c)
                _lepton_pdgid.append(13)
        else:
            pass
        if len(lepton.electron_pt[i]) == 0 and len(lepton.Muon_pt[i]) == 0:
            _lepton_pt.append(0)
            _lepton_eta.append(0)
            _lepton_phi.append(0)
            _lepton_pdgid.append(0)
        else: pass
        tmp_lepton_pt.append(np.asanyarray(_lepton_pt, dtype=object))
        tmp_lepton_eta.append(np.asanyarray(_lepton_eta, dtype=object))
        tmp_lepton_phi.append(np.asanyarray(_lepton_phi, dtype=object))
        tmp_lepton_pdgid.append(np.asanyarray(_lepton_pdgid, dtype=object))

    #Mark which jet in each event pass the selection.
    print("+------------------------------------------------------------------------------------------------------+")
    print("Start jet selection.")
    print("+------------------------------------------------------------------------------------------------------+")
    marker_event, marker_jet, marker_btag = event_selection(jet.pt, jet.eta, jet.phi, jet.btag, tmp_lepton_pt, tmp_lepton_eta, tmp_lepton_phi, MODEL)
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

    if MODEL == 'ttbar_lep_left' or MODEL == "ttbar_lep_right":
        lepton_pt = []
        lepton_eta = []
        lepton_phi = []
        lepton_pdgid = []
        for i in tqdm.trange(len(tmp_lepton_pt)):
            if marker_event[i] == 1:
                if MODEL == 'ttbar_lep_left' or MODEL == "ttbar_lep_right":
                    lepton_pt_tmp = []
                    lepton_eta_tmp = []
                    lepton_phi_tmp = []
                    lepton_pdgid_tmp = []
                    for j in range(len(tmp_lepton_pt[i])):
                        lepton_pt_tmp.append(tmp_lepton_pt[i][j])
                        lepton_eta_tmp.append(tmp_lepton_eta[i][j])
                        lepton_phi_tmp.append(tmp_lepton_phi[i][j])
                        lepton_pdgid_tmp.append(tmp_lepton_pdgid[i][j])
                    lepton_pt.append(np.asanyarray(lepton_pt_tmp, dtype=object))
                    lepton_eta.append(np.asanyarray(lepton_eta_tmp, dtype=object))
                    lepton_phi.append(np.asanyarray(lepton_phi_tmp, dtype=object))
                    lepton_pdgid.append(np.asanyarray(lepton_pdgid_tmp, dtype=object))
        lepton_pt = np.asanyarray(lepton_pt, dtype=object)
        lepton_eta = np.asanyarray(lepton_eta, dtype=object)
        lepton_phi = np.asanyarray(lepton_phi, dtype=object)
        lepton_pdgid = np.asanyarray(lepton_pdgid, dtype=object)

    for i in tqdm.trange(len(jet.event)):
        if marker_event[i] == 1:
            jet_pt_tmp = []
            jet_eta_tmp = []
            jet_phi_tmp = []
            jet_mass_tmp = []
            jet_btag_tmp = []
            
            if MODEL == 'ttbar_lep_left' or MODEL == "ttbar_lep_right":
                #Lepton signal removal for ttbar semi-leptonic model
                if len(tmp_lepton_eta[i]) != 1:
                    print(i, len(tmp_lepton_eta[i]), len(tmp_lepton_eta[i]), len(tmp_lepton_phi[i]), tmp_lepton_pt[i], tmp_lepton_eta[i])
                for j in range(len(jet.pt[i])):
                    if delta_R(jet.eta[i][j], jet.phi[i][j], tmp_lepton_eta[i][0], tmp_lepton_phi[i][0]) > 0.4:
                        jet_pt_tmp.append(jet.pt[i][j])
                        jet_eta_tmp.append(jet.eta[i][j])
                        jet_phi_tmp.append(jet.phi[i][j])
                        jet_mass_tmp.append(jet.mass[i][j])
                        jet_btag_tmp.append(jet.btag[i][j])
                    else : pass 
            else: 
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
        if MODEL == 'ttbar_lep_left' or MODEL == "ttbar_lep_right":
            lepton_pt = np.asanyarray(lepton_pt, dtype=object)
            lepton_eta = np.asanyarray(lepton_eta, dtype=object)
            lepton_phi = np.asanyarray(lepton_phi, dtype=object)
            lepton_pdgid = np.asanyarray(lepton_pdgid, dtype=object)

    MET = []
    MET_ETA = []
    MET_PHI = []
    for i in range(len(missing_et.MET)):
        if marker_event[i] == 1:
             for a, b, c in zip(missing_et.MET[i], missing_et.eta[i], missing_et.phi[i]):
                MET.append(a)
                MET_ETA.append(b)
                MET_PHI.append(c)
    MET = np.asanyarray(MET, dtype=object)
    MET_ETA = np.asanyarray(MET_ETA, dtype=object)
    MET_PHI = np.asanyarray(MET_PHI, dtype=object)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Finished to record the kinematics variables of jets in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    
    print("+------------------------------------------------------------------------------------------------------+")
    print("Starting parton tracing and looking for its daughter.")
    print("+------------------------------------------------------------------------------------------------------+")
    #Particle tracing and daughter finding section
    start = time.time()
    if MODEL == 'ttbar' or MODEL == 'ttbar_lep_left' or MODEL =='ttbar_lep_right':

        top_idx, top_daughter_idx_1, top_daughter_pid_1, top_daughter_idx_2, top_daughter_pid_2 = [], [], [], [], []
        top_bar_idx, top_bar_daughter_idx_1, top_bar_daughter_pid_1, top_bar_daughter_idx_2, top_bar_daughter_pid_2 = [], [], [], [], []
        _src_top, _src_anti_top, _index  = [], [], []
        for i in range(len(particle.event)):
            if marker_event[i] == 1:
                _index.append(i)
                _src_top.append([particle.dataframelize(i), PID.top, 22, MODEL])
                _src_anti_top.append([particle.dataframelize(i), PID.anti_top, 22, MODEL])
        print("Using {0} process for accelerating speed.".format(PROCESS))
        with mp.Pool(PROCESS) as p:
            _result_top = p.starmap(particle_tracing, _src_top)
            p.close()
            p.join()
        print("Top tracing finished.")
        with mp.Pool(PROCESS) as p:
            _result_anti_top = p.starmap(particle_tracing, _src_anti_top)
            p.close()
            p.join()
        print("Anti-Top tracing finished.")
        for i in range(len(_result_top)):
            top_idx.append(_result_top[i][0])
            top_daughter_idx_1.append(_result_top[i][1])
            top_daughter_pid_1.append(_result_top[i][2])
            top_daughter_idx_2.append(_result_top[i][3])
            top_daughter_pid_2.append(_result_top[i][4])
            top_bar_idx.append(_result_anti_top[i][0])
            top_bar_daughter_idx_1.append(_result_anti_top[i][1])
            top_bar_daughter_pid_1.append(_result_anti_top[i][2])
            top_bar_daughter_idx_2.append(_result_anti_top[i][3])
            top_bar_daughter_pid_2.append(_result_anti_top[i][4])

        _src_top_d, _src_anti_top_d = [], []
        parton_array = np.zeros([ len(_index) , NUM_OF_DAUGHTER, 7])

        for i in range(len(_index)):
            j = _index[i]
            _src_top_d.append([particle.dataframelize(j), top_daughter_idx_1[i], top_daughter_idx_2[i]])
            _src_anti_top_d.append([particle.dataframelize(j), top_bar_daughter_idx_1[i], top_bar_daughter_idx_2[i]])
        with mp.Pool(PROCESS) as p:
            _result_top = p.starmap(quark_finder, _src_top_d)
            p.close()
            p.join()
        print("Daughter of Top's daughter found.")
        with mp.Pool(PROCESS) as p:
            _result_anti_top = p.starmap(quark_finder, _src_anti_top_d)
            p.close()
            p.join()
        print("Daughter of Anti-Top's daughter found.")
        
        for i in range(len(_index)):
            parton_array[i][0][0], parton_array[i][1][0], parton_array[i][2][0] = _result_top[i][0], _result_top[i][1], _result_top[i][2]
            parton_array[i][3][0], parton_array[i][4][0], parton_array[i][5][0] = _result_anti_top[i][0], _result_anti_top[i][1], _result_anti_top[i][2]
        print("+------------------------------------------------------------------------------------------------------+")
        print("Parton tracing section complete. The daughter of W+/W- and bbbar has been found. Cost: {0:.1f} s".format(time.time()-start))
        print("+------------------------------------------------------------------------------------------------------+")
    elif MODEL == 'ttH':
        top_idx, top_daughter_idx_1, top_daughter_pid_1, top_daughter_idx_2, top_daughter_pid_2 = [], [], [], [], []
        top_bar_idx, top_bar_daughter_idx_1, top_bar_daughter_pid_1, top_bar_daughter_idx_2, top_bar_daughter_pid_2 = [], [], [], [], []
        higgs_idx, higgs_daughter_idx_1, higgs_daughter_pid_1, higgs_daughter_idx_2, higgs_daughter_pid_2 = [], [], [], [], []
        _src_top, _src_anti_top, _src_higgs, _index  = [], [], [], []
        
        for i in range(len(particle.event)):
            if marker_event[i] == 1:
                _index.append(i)
                _src_top.append([particle.dataframelize(i), PID.top, 22, MODEL])
                _src_anti_top.append([particle.dataframelize(i), PID.anti_top, 22, MODEL])
                _src_higgs.append([particle.dataframelize(i), PID.higgs, 22, MODEL])
        print("Using {0} process for accelerating speed.".format(PROCESS))
        with mp.Pool(PROCESS) as p:
            _result_top = p.starmap(particle_tracing, _src_top)
            p.close()
            p.join()
        print("Top tracing finished.")
        with mp.Pool(PROCESS) as p:
            _result_anti_top = p.starmap(particle_tracing, _src_anti_top)
            p.close()
            p.join()
        print("Anti-Top tracing finished.")
        with mp.Pool(PROCESS) as p:
            _result_h = p.starmap(particle_tracing, _src_higgs)
            p.close()
            p.join()
        print("Higgs tracing finished.")
        
        for i in range(len(_index)):
            top_idx.append(_result_top[i][0])
            top_daughter_idx_1.append(_result_top[i][1])
            top_daughter_pid_1.append(_result_top[i][2])
            top_daughter_idx_2.append(_result_top[i][3])
            top_daughter_pid_2.append(_result_top[i][4])
            top_bar_idx.append(_result_anti_top[i][0])
            top_bar_daughter_idx_1.append(_result_anti_top[i][1])
            top_bar_daughter_pid_1.append(_result_anti_top[i][2])
            top_bar_daughter_idx_2.append(_result_anti_top[i][3])
            top_bar_daughter_pid_2.append(_result_anti_top[i][4])
            higgs_idx.append(_result_h[i][0])
            higgs_daughter_idx_1.append(_result_h[i][1])
            higgs_daughter_pid_1.append(_result_h[i][2])
            higgs_daughter_idx_2.append(_result_h[i][3])
            higgs_daughter_pid_2.append(_result_h[i][4])
        
        _src_top_d, _src_anti_top_d  = [], []
        parton_array = np.zeros([ len(_index) , NUM_OF_DAUGHTER, 7])

        for i in range(len(_index)):
            j = _index[i]
            _src_top_d.append([particle.dataframelize(j), top_daughter_idx_1[i], top_daughter_idx_2[i]])
            _src_anti_top_d.append([particle.dataframelize(j), top_bar_daughter_idx_1[i], top_bar_daughter_idx_2[i]])
        with mp.Pool(PROCESS) as p:
            _result_top = p.starmap(quark_finder, _src_top_d)
            p.close()
            p.join()
        print("Daughter of Top's daughter found.")
        with mp.Pool(PROCESS) as p:
            _result_anti_top = p.starmap(quark_finder, _src_anti_top_d)
            p.close()
            p.join()
        print("Daughter of Anti-Top's daughter found.")
        for i in range(len(_index)):
            parton_array[i][0][0], parton_array[i][1][0], parton_array[i][2][0] = _result_top[i][0], _result_top[i][1], _result_top[i][2]
            parton_array[i][3][0], parton_array[i][4][0], parton_array[i][5][0] = _result_anti_top[i][0], _result_anti_top[i][1], _result_anti_top[i][2]
            parton_array[i][6][0], parton_array[i][7][0] = higgs_daughter_idx_1[i], higgs_daughter_idx_2[i]
        print("+------------------------------------------------------------------------------------------------------+")
        print("Parton tracing section complete. The daughter of W+/W-, bbbar, and Higgs has been found. Cost: {0:.1f} s".format(time.time()-start))
        print("+------------------------------------------------------------------------------------------------------+")
    elif MODEL == 'four_top':
        top_1_idx, top_1_daughter_idx_1, top_1_daughter_pid_1, top_1_daughter_idx_2, top_1_daughter_pid_2 = [], [], [], [], []
        top_2_idx, top_2_daughter_idx_1, top_2_daughter_pid_1, top_2_daughter_idx_2, top_2_daughter_pid_2 = [], [], [], [], []
        top_1_bar_idx, top_1_bar_daughter_idx_1, top_1_bar_daughter_pid_1, top_1_bar_daughter_idx_2, top_1_bar_daughter_pid_2 = [], [], [], [], []
        top_2_bar_idx, top_2_bar_daughter_idx_1, top_2_bar_daughter_pid_1, top_2_bar_daughter_idx_2, top_2_bar_daughter_pid_2 = [], [], [], [], []
        
        _src_top, _src_anti_top, _index  = [], [], []
        for i in range(len(particle.event)):
            if marker_event[i] == 1:
                _index.append(i)
                _src_top.append([particle.dataframelize(i), PID.top, 22, MODEL])
                _src_anti_top.append([particle.dataframelize(i), PID.anti_top, 22, MODEL])
        print("Using {0} process for accelerating speed.".format(PROCESS))
        with mp.Pool(PROCESS) as p:
            _result_top = p.starmap(particle_tracing, _src_top)
            p.close()
            p.join()
        print("Top tracing finished.")
        with mp.Pool(PROCESS) as p:
            _result_anti_top = p.starmap(particle_tracing, _src_anti_top)
            p.close()
            p.join()
        print("Anti-Top tracing finished.")
        for i in range(len(_index)):
            top_1_idx.append(_result_top[i][0])
            top_2_idx.append(_result_top[i][1])
            top_1_daughter_idx_1.append(_result_top[i][2])
            top_1_daughter_pid_1.append(_result_top[i][3]) 
            top_1_daughter_idx_2.append(_result_top[i][4])
            top_1_daughter_pid_2.append(_result_top[i][5])
            top_2_daughter_idx_1.append(_result_top[i][6]) 
            top_2_daughter_pid_1.append(_result_top[i][7])
            top_2_daughter_idx_2.append(_result_top[i][8])
            top_2_daughter_pid_2.append(_result_top[i][9])
            top_1_bar_idx.append(_result_anti_top[i][0])
            top_2_bar_idx.append(_result_anti_top[i][1])
            top_1_bar_daughter_idx_1.append(_result_anti_top[i][2]) 
            top_1_bar_daughter_pid_1.append(_result_anti_top[i][3])
            top_1_bar_daughter_idx_2.append(_result_anti_top[i][4])
            top_1_bar_daughter_pid_2.append(_result_anti_top[i][5])
            top_2_bar_daughter_idx_1.append(_result_anti_top[i][6]) 
            top_2_bar_daughter_pid_1.append(_result_anti_top[i][7])
            top_2_bar_daughter_idx_2.append(_result_anti_top[i][8])
            top_2_bar_daughter_pid_2.append(_result_anti_top[i][9])
        
        _src_top_d_1, _src_top_d_2, _src_anti_top_d_1, _src_anti_top_d_2 = [], [], [], []
        
        parton_array = np.zeros([ len(_index) , NUM_OF_DAUGHTER, 7])

        for i in range(len(_index)):
            j = _index[i]
            _src_top_d_1.append([particle.dataframelize(j), top_1_daughter_idx_1[i], top_1_daughter_idx_2[i]])
            _src_top_d_2.append([particle.dataframelize(j), top_2_daughter_idx_1[i], top_2_daughter_idx_2[i]])
            _src_anti_top_d_1.append([particle.dataframelize(j), top_1_bar_daughter_idx_1[i], top_1_bar_daughter_idx_2[i]])
            _src_anti_top_d_2.append([particle.dataframelize(j), top_2_bar_daughter_idx_1[i], top_2_bar_daughter_idx_2[i]])
        with mp.Pool(PROCESS) as p:
            _result_top_1 = p.starmap(quark_finder, _src_top_d_1)
            p.close()
            p.join()
        print("Daughter of Top_1's daughter found.") 
        with mp.Pool(PROCESS) as p:
            _result_top_2 = p.starmap(quark_finder, _src_top_d_2)
            p.close()
            p.join()
        print("Daughter of Top_2's daughter found.") 
        with mp.Pool(PROCESS) as p:
            _result_anti_top_1 = p.starmap(quark_finder, _src_anti_top_d_1)
            p.close()
            p.join()
        print("Daughter of Anti-Top_1's daughter found.") 
        with mp.Pool(PROCESS) as p:
            _result_anti_top_2 = p.starmap(quark_finder, _src_anti_top_d_2)
            p.close()
            p.join()
        print("Daughter of Anti-Top_2's daughter found.") 
        for i in range(len(_index)):
            parton_array[i][0][0], parton_array[i][1][0], parton_array[i][2][0] = _result_top_1[i][0], _result_top_1[i][1], _result_top_1[i][2]
            parton_array[i][3][0], parton_array[i][4][0], parton_array[i][5][0] = _result_top_2[i][0], _result_top_2[i][1], _result_top_2[i][2]
            parton_array[i][6][0], parton_array[i][7][0], parton_array[i][8][0] = _result_anti_top_1[i][0], _result_anti_top_1[i][1], _result_anti_top_1[i][2]
            parton_array[i][9][0], parton_array[i][10][0], parton_array[i][11][0] = _result_anti_top_2[i][0], _result_anti_top_2[i][1], _result_anti_top_2[i][2]
        print("+------------------------------------------------------------------------------------------------------+")
        print("Parton tracing section complete. The daughter of W+/W- and bbbar has been found. Cost: {0:.1f} s".format(time.time()-start))
        print("+------------------------------------------------------------------------------------------------------+")
    else :
        print("Please select a correct model.")


    

    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording the kinematics variables of partons in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    parton_pdgid = []
    parton_barcode = []
    parton_pt = []
    parton_eta = []
    parton_phi = []
    parton_mass = []


    for i in tqdm.trange(len(_index)):
        _parton_pdgid = []
        _parton_barcode = []
        _parton_pt = []
        _parton_eta = []
        _parton_phi = []
        _parton_mass = []
        for j in range(NUM_OF_DAUGHTER):
            k = _index[i]
            dataset = particle.dataframelize(k)

            _parton_pdgid.append(dataset.iloc[int(parton_array[i][j][0]), 6])
            _parton_barcode.append(barcode[j])
            _parton_pt.append(dataset.iloc[int(parton_array[i][j][0]), 7])
            _parton_eta.append(dataset.iloc[int(parton_array[i][j][0]), 8])
            _parton_phi.append(dataset.iloc[int(parton_array[i][j][0]), 9])
            _parton_mass.append(dataset.iloc[int(parton_array[i][j][0]), 10])


        parton_pdgid.append(np.asanyarray(_parton_pdgid, dtype=object))
        parton_barcode.append(np.asanyarray(_parton_barcode, dtype=object))
        parton_pt.append(np.asanyarray(_parton_pt, dtype=object))
        parton_eta.append(np.asanyarray(_parton_eta, dtype=object))
        parton_phi.append(np.asanyarray(_parton_phi, dtype=object))
        parton_mass.append(np.asanyarray(_parton_mass, dtype=object))
    parton_pdgid = np.asanyarray(parton_pdgid, dtype=object)
    parton_barcode = np.asanyarray(parton_barcode, dtype=object)
    parton_pt = np.asanyarray(parton_pt, dtype=object)
    parton_eta = np.asanyarray(parton_eta, dtype=object)
    parton_phi = np.asanyarray(parton_phi, dtype=object)
    parton_mass = np.asanyarray(parton_mass, dtype=object)
    
    if MODEL == 'ttbar_lep_left':
        print("Recording simulation lepton kinematic properties.")
        simulation_lepton_pdgid = np.zeros(len(_index))
        simulation_lepton_barcode = np.zeros(len(_index))
        simulation_lepton_pt = np.zeros(len(_index))
        simulation_lepton_eta = np.zeros(len(_index))
        simulation_lepton_phi = np.zeros(len(_index))
        simulation_lepton_mass = np.zeros(len(_index))
        simulation_neutrino_pdgid = np.zeros(len(_index))
        simulation_neutrino_barcode = np.zeros(len(_index))
        simulation_neutrino_pt = np.zeros(len(_index))
        simulation_neutrino_eta = np.zeros(len(_index))
        simulation_neutrino_phi = np.zeros(len(_index))
        simulation_neutrino_mass = np.zeros(len(_index))
        for i in tqdm.trange(len(_index)):
            for j in range(1,3):
                if parton_pdgid[i][j] == -11 or parton_pdgid[i][j] == -13:
                    simulation_lepton_pdgid[i] = parton_pdgid[i][j]
                    simulation_lepton_barcode[i] = parton_barcode[i][j]
                    simulation_lepton_pt[i] = parton_pt[i][j]
                    simulation_lepton_eta[i] = parton_eta[i][j]
                    simulation_lepton_phi[i] = parton_phi[i][j]
                    simulation_lepton_mass[i] = parton_mass[i][j]
                    
                else: 
                    simulation_neutrino_pdgid[i] = parton_pdgid[i][j]
                    simulation_neutrino_barcode[i] = parton_barcode[i][j]
                    simulation_neutrino_pt[i] = parton_pt[i][j]
                    simulation_neutrino_eta[i] = parton_eta[i][j]
                    simulation_neutrino_phi[i] = parton_phi[i][j]
                    simulation_neutrino_mass[i] = parton_mass[i][j]

        parton_pdgid = np.delete(parton_pdgid, [1,2], 1)
        parton_barcode = np.delete(parton_barcode, [1,2], 1)
        parton_pt = np.delete(parton_pt, [1,2], 1)
        parton_eta = np.delete(parton_eta, [1,2], 1)
        parton_phi = np.delete(parton_phi, [1,2], 1)
        parton_mass = np.delete(parton_mass, [1,2], 1)

    if MODEL == "ttbar_lep_right":
        print("Recording simulation lepton kinematic properties.")
        simulation_lepton_pdgid = np.zeros(len(_index))
        simulation_lepton_barcode = np.zeros(len(_index))
        simulation_lepton_pt = np.zeros(len(_index))
        simulation_lepton_eta = np.zeros(len(_index))
        simulation_lepton_phi = np.zeros(len(_index))
        simulation_lepton_mass = np.zeros(len(_index))
        simulation_neutrino_pdgid = np.zeros(len(_index))
        simulation_neutrino_barcode = np.zeros(len(_index))
        simulation_neutrino_pt = np.zeros(len(_index))
        simulation_neutrino_eta = np.zeros(len(_index))
        simulation_neutrino_phi = np.zeros(len(_index))
        simulation_neutrino_mass = np.zeros(len(_index))
        for i in tqdm.trange(len(_index)):
            for j in range(4,6):
                if parton_pdgid[i][j] == 11 or parton_pdgid[i][j] == 13:
                    simulation_lepton_pdgid[i] = parton_pdgid[i][j]
                    simulation_lepton_barcode[i] = parton_barcode[i][j]
                    simulation_lepton_pt[i] = parton_pt[i][j]
                    simulation_lepton_eta[i] = parton_eta[i][j]
                    simulation_lepton_phi[i] = parton_phi[i][j]
                    simulation_lepton_mass[i] = parton_mass[i][j]
                else: 
                    simulation_neutrino_pdgid[i] = parton_pdgid[i][j]
                    simulation_neutrino_barcode[i] = parton_barcode[i][j]
                    simulation_neutrino_pt[i] = parton_pt[i][j]
                    simulation_neutrino_eta[i] = parton_eta[i][j]
                    simulation_neutrino_phi[i] = parton_phi[i][j]
                    simulation_neutrino_mass[i] = parton_mass[i][j]
        parton_pdgid = np.delete(parton_pdgid, [4, 5], 1)
        parton_barcode = np.delete(parton_barcode, [4, 5], 1)
        parton_pt = np.delete(parton_pt, [4, 5], 1)
        parton_eta = np.delete(parton_eta, [4, 5], 1)
        parton_phi = np.delete(parton_phi, [4, 5], 1)
        parton_mass = np.delete(parton_mass, [4, 5], 1)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Finished to record the kinematics variables of partons in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")

    
    print("+------------------------------------------------------------------------------------------------------+")
    print("Starting parton-jet matching.")
    print("+------------------------------------------------------------------------------------------------------+")
    jet_parton_index = []
    parton_jet_index = []
    _src_delta_R = []
    start = time.time()
    for i in tqdm.trange(len(jet_pt)):
        _src_delta_R.append([NUM_OF_PARTON, len(jet_pt[i]), parton_eta[i], parton_phi[i], jet_eta[i], jet_phi[i], 0.4, MODEL])
    print("Using {0} process for accelerating speed.".format(PROCESS))
    with mp.Pool(PROCESS) as p:
        _result_delta_R = p.starmap(deltaR_matching, _src_delta_R)
        p.close()
        p.join()
    
    for i in range(len(_result_delta_R)):
        jet_parton_index.append(_result_delta_R[i][0])
        parton_jet_index.append(_result_delta_R[i][1])
    print("+------------------------------------------------------------------------------------------------------+")
    print("Parton-jet matching finished. Cost: {0:.1f} s".format(time.time()-start))
    print("+------------------------------------------------------------------------------------------------------+")
    jet_parton_index = np.asanyarray(jet_parton_index, dtype=object)
    parton_jet_index = np.asanyarray(parton_jet_index, dtype=object)
    
    if MODEL == 'ttbar_lep_left' or MODEL == "ttbar_lep_right":
        print("+------------------------------------------------------------------------------------------------------+")
        print("Starting lepton matching.")
        print("+------------------------------------------------------------------------------------------------------+")
        lepton_delta_R_result = np.zeros(len(simulation_lepton_pt))
        for i in range(len(simulation_lepton_pt)):
            _delta_R = delta_R(simulation_lepton_eta[i], simulation_lepton_phi[i], lepton_eta[i][0], lepton_phi[i][0])
            if _delta_R < 0.4:
                lepton_delta_R_result[i] = 1
            else : 
                lepton_delta_R_result[i] = 0
        print("+------------------------------------------------------------------------------------------------------+")
        print("Lepton matching finished.")
        print("+------------------------------------------------------------------------------------------------------+") 
        
        print("+------------------------------------------------------------------------------------------------------+")
        print("Starting neutrino matching.")
        print("+------------------------------------------------------------------------------------------------------+")
        MET_delta_R_result = np.zeros(len(simulation_neutrino_eta))
        for i in range(len(simulation_neutrino_eta)):
            _delta_R = deltaPhi(simulation_neutrino_phi[i], MET_PHI[i])
            if np.abs(_delta_R) < 0.4:
                MET_delta_R_result[i] = 1
            else : 
                MET_delta_R_result[i] = 0
        print("+------------------------------------------------------------------------------------------------------+")
        print("Neutrino matching finished.")
        print("+------------------------------------------------------------------------------------------------------+") 

    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording barcode information.")
    print("+------------------------------------------------------------------------------------------------------+")

    jet_barcode = []
    for i in tqdm.trange(len(jet_pt)):
        _jet_barcode = barcode_recorder(jet_parton_index[i], MODEL)
        jet_barcode.append(_jet_barcode)
    jet_barcode = np.asanyarray(jet_barcode, dtype=object)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Barcode information has beed record.")
    print("+------------------------------------------------------------------------------------------------------+")

    if MODEL == 'ttH':
        target = [i for i in range(6)]
        N_match_top_in_event = np.zeros([len(jet_pt)])
        for i in tqdm.trange(len(jet_parton_index)):
            intersetion = set(target).intersection(jet_parton_index[i])
            if intersetion.intersection(set([0, 1, 2])) == {0,1,2} and intersetion.intersection(set([3, 4, 5])) == {3,4,5} :
                N_match_top_in_event[i] = 2
            elif intersetion.intersection(set([0, 1, 2])) != {0,1,2} or intersetion.intersection(set([3, 4, 5])) != {3,4,5}:
                N_match_top_in_event[i] = 1
            elif intersetion.intersection(set([0, 1, 2])) != {0,1,2} and intersetion.intersection(set([3, 4, 5])) != {3,4,5}:
                N_match_top_in_event[i] = 0
            else : pass

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
            if intersetion.intersection(set([0, 1, 2])) == {0,1,2} and intersetion.intersection(set([3, 4, 5])) == {3,4,5} :
                N_match_top_in_event[i] = 2
            elif intersetion.intersection(set([0, 1, 2])) != {0,1,2} or intersetion.intersection(set([3, 4, 5])) != {3,4,5}:
                N_match_top_in_event[i] = 1
            elif intersetion.intersection(set([0, 1, 2])) != {0,1,2} and intersetion.intersection(set([3, 4, 5])) != {3,4,5}:
                N_match_top_in_event[i] = 0
        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2)  ))
        print("+------------------------------------------------------------------------------------------------------+")
    elif MODEL == 'four_top':

        target = [i for i in range(NUM_OF_PARTON)]
        N_match_top_in_event = np.zeros([len(jet_pt)])
        for i in tqdm.trange(len(jet_parton_index)):
            
            intersetion = set(target).intersection(jet_parton_index[i])
            count_inter = 0
            if intersetion.intersection(set([0, 1, 2])) == {0,1,2}:
                count_inter += 1
            if intersetion.intersection(set([3, 4, 5])) == {3,4,5}:
                count_inter += 1
            if intersetion.intersection(set([6, 7, 8])) == {6,7,8}:
                count_inter += 1
            if intersetion.intersection(set([9, 10, 11])) == {9,10,11}:
                count_inter += 1

            N_match_top_in_event[i] = count_inter
        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.\nFound {2} events with 3 ttbar candidate exist.\nFound {3} events with 4 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2), np.sum(N_match_top_in_event == 3), np.sum(N_match_top_in_event == 4)  ))
        print("+------------------------------------------------------------------------------------------------------+")
    elif MODEL == 'ttbar_lep_left' or MODEL == 'ttbar_lep_right':
        N_match_top_in_event = np.zeros([len(jet_pt)])
        target = [i for i in range(NUM_OF_PARTON)]
        for i in tqdm.trange(len(jet_parton_index)):
            intersetion = set(target).intersection(jet_parton_index[i])
            if MET_delta_R_result[i] == 1 and lepton_delta_R_result[i] == 1:              
                if len(intersetion) == 4:
                    N_match_top_in_event[i] = 2
                elif 3 in intersetion:
                    N_match_top_in_event[i] = 1
                elif intersetion.intersection(set([0, 1, 2])) == {0,1,2} and (3 in intersetion) == False:
                    N_match_top_in_event[i] = 1
            else: 
                if intersetion.intersection(set([0, 1, 2])) == {0, 1, 2}:
                    N_match_top_in_event[i] = 1
                else: pass
        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2) ))
        print("+------------------------------------------------------------------------------------------------------+")
    else : pass
    

    print("+------------------------------------------------------------------------------------------------------+")
    print("Writing event record to the npz file.")
    print("+------------------------------------------------------------------------------------------------------+")

    if MODEL == 'ttbar' or MODEL == 'four_top':
        np.savez_compressed(OUTPUT_FILE, 
                            jet_parton_index=jet_parton_index,
                            jet_barcode=jet_barcode,
                            jet_pt=jet_pt,
                            jet_eta=jet_eta,
                            jet_phi=jet_phi,
                            jet_mass=jet_mass,
                            jet_btag=jet_btag,
                            parton_jet_index=parton_jet_index,
                            parton_pdgid=parton_pdgid,
                            parton_barcode=parton_barcode,
                            parton_pt=parton_pt,
                            parton_eta=parton_eta,
                            parton_phi=parton_phi,
                            parton_mass=parton_mass,
                            N_match_top_in_event=N_match_top_in_event)
    elif MODEL == 'ttH':
        np.savez_compressed(OUTPUT_FILE, 
                            jet_parton_index=jet_parton_index,
                            jet_barcode=jet_barcode,
                            jet_pt=jet_pt,
                            jet_eta=jet_eta,
                            jet_phi=jet_phi,
                            jet_mass=jet_mass,
                            jet_btag=jet_btag,
                            parton_jet_index=parton_jet_index,
                            parton_pdgid=parton_pdgid,
                            parton_barcode=parton_barcode,
                            parton_pt=parton_pt,
                            parton_eta=parton_eta,
                            parton_phi=parton_phi,
                            parton_mass=parton_mass,
                            N_match_top_in_event=N_match_top_in_event,
                            N_match_higgs_in_event=N_match_higgs_in_event)
    elif MODEL == 'ttbar_lep_left' or MODEL == 'ttbar_lep_right':
        np.savez_compressed(OUTPUT_FILE, 
                            jet_parton_index=jet_parton_index,
                            jet_barcode=jet_barcode,
                            jet_pt=jet_pt,
                            jet_eta=jet_eta,
                            jet_phi=jet_phi,
                            jet_mass=jet_mass,
                            jet_btag=jet_btag,
                            parton_jet_index=parton_jet_index,
                            parton_pdgid=parton_pdgid,
                            parton_barcode=parton_barcode,
                            parton_pt=parton_pt,
                            parton_eta=parton_eta,
                            parton_phi=parton_phi,
                            parton_mass=parton_mass,
                            N_match_top_in_event=N_match_top_in_event,
                            lepton_pt=lepton_pt,
                            lepton_eta=lepton_eta,
                            lepton_phi=lepton_phi,
                            lepton_pdgid=lepton_pdgid,
                            MET=MET,
                            MET_ETA=MET_ETA,
                            MET_PHI=MET_PHI,
                            simulation_neutrino_pt=simulation_neutrino_pt,
                            simulation_neutrino_eta=simulation_neutrino_eta,
                            simulation_neutrino_phi=simulation_neutrino_phi,
                            simulation_neutrino_pdgid=simulation_neutrino_pdgid,
                            simulation_neutrino_barcode=simulation_neutrino_barcode,
                            simulation_neutrino_mass=simulation_neutrino_mass,
                            simulation_lepton_pt=simulation_lepton_pt,
                            simulation_lepton_eta=simulation_lepton_eta,
                            simulation_lepton_phi=simulation_lepton_phi,
                            simulation_lepton_pdgid=simulation_lepton_pdgid,
                            simulation_lepton_mass=simulation_lepton_mass,
                            simulation_lepton_barcode=simulation_lepton_barcode)

    print("+------------------------------------------------------------------------------------------------------+")
    print("Event record has been send to {0}.".format(OUTPUT_FILE))
    print("+------------------------------------------------------------------------------------------------------+")
