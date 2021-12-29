"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan
Mail: davidho@gapp.nthu.edu.tw
Modified by Hideki Okawa (Fudan U, hideki.okawa@cern.ch) and Mike Fenton (UCI, m.fenton@uci.edu)
"""
import numpy as np
import itertools, uproot, sys, os, tqdm
import pandas as pd
from typing import Union
import numba as nb
import h5py
from collections import OrderedDict

class pdgid():
    """
    This is a function for configurating PDG id.
    """
    def __init__(self):
        self.w_plus = 24
        self.w_minus = -24
        self.down = 1
        self.anti_down = -1 
        self.up = 2
        self.anti_up = -2 
        self.strange = 3 
        self.anti_strange = -3
        self.charm = 4
        self.anti_charm = -4
        self.bottom = 5 
        self.anti_bottom = -5
        self.top = 6
        self.anti_top = -6
        self.higgs = 25
        self.electron = 11 
        self.positron = -11
        self.electron_neutrino = 12
        self.anti_electron_neutrino = -12
        self.muon = 13
        self.anti_muon = -13
        self.muon_neutrino = 14
        self.anti_muon_neutrino = -14
        self.tau = 17 
        self.anti_tau = -17
        self.tau_neutrino = 18
        self.anti_tau_neutrino = -18
        self.z_plus = 23
        self.z_minus = -23
        self.w_plus = 24
        self.w_minus = -24

PID = pdgid()

class IO_module:
    """
    This is a I/O module for HEP-jet-assignment project.
    This module can help to load the data from root files/npz files and output as npz/hdf5 files.
    """
    def __init__(self, PATH, MODEL, MULTI=False, **kargs):
        self.path = PATH
        self.multi = MULTI
        self.model = MODEL
        self.require_lepton = ["ttbar_lep", "ttbar_lep_left", "ttbar_lep_right"]
        self.kargs = kargs
    def read_ROOT(self) -> dict:
        # If loading multi-root files, using this function to concatenate dataset.
        if self.multi:
            files_under_path = os.listdir(self.path)
            file_list = [os.path.join(self.path, a) for a in files_under_path]
            num_of_files = len(file_list)
            _data = []
            count = 1
            for a in file_list:
                try:
                    print(f"Padding root file from {a}. Progress: {count}/{num_of_files}.")
                    tmp = uproot.open(a)['Delphes']
                    _data.append(tmp)
                    count += 1
                except:
                    print('Please check input file path.')
            for i in tqdm.trange(num_of_files):
                if i == 0 :
                    _particle_event = _data[i].array('Event')
                    _particle_pt = _data[i].array('Particle.PT')
                    _particle_eta = _data[i].array('Particle.Eta')
                    _particle_phi = _data[i].array('Particle.Phi')
                    _particle_pid = _data[i].array('Particle.PID')
                    _particle_M1 = _data[i].array('Particle.M1')
                    _particle_M2 = _data[i].array('Particle.M2')
                    _particle_D1 = _data[i].array('Particle.D1')
                    _particle_D2 = _data[i].array('Particle.D2')
                    _particle_status = _data[i].array('Particle.Status')
                    _particle_rapidity = _data[i].array('Particle.Rapidity')
                    _particle_mass = _data[i].array('Particle.Mass')
                    _particle_charge = _data[i].array('Particle.Charge')

                    _jet_event = _data[i].array('Event')
                    _jet_pt = _data[i].array('Jet.PT')
                    _jet_eta = _data[i].array('Jet.Eta')
                    _jet_phi = _data[i].array('Jet.Phi')
                    _jet_btag = _data[i].array('Jet.BTag')
                    _jet_mass = _data[i].array('Jet.Mass')
                    _jet_charge = _data[i].array('Jet.Charge')
                    _num_of_jets = _data[i].array('Jet')
                else: 
                    _particle_event = np.concatenate((_particle_event, _data[i].array('Event')))
                    _particle_pt = np.concatenate((_particle_pt, _data[i].array('Particle.PT')))
                    _particle_eta = np.concatenate((_particle_eta, _data[i].array('Particle.Eta')))
                    _particle_phi = np.concatenate((_particle_phi, _data[i].array('Particle.Phi')))
                    _particle_pid = np.concatenate((_particle_pid,_data[i].array('Particle.PID')))
                    _particle_M1 = np.concatenate((_particle_M1, _data[i].array('Particle.M1')))
                    _particle_M2 = np.concatenate((_particle_M2, _data[i].array('Particle.M2')))
                    _particle_D1 = np.concatenate((_particle_D1, _data[i].array('Particle.D1')))
                    _particle_D2 = np.concatenate((_particle_D2, _data[i].array('Particle.D2')))
                    _particle_status = np.concatenate((_particle_status, _data[i].array('Particle.Status')))
                    _particle_rapidity = np.concatenate((_particle_rapidity, _data[i].array('Particle.Rapidity')))
                    _particle_mass = np.concatenate((_particle_mass, _data[i].array('Particle.Mass')))
                    _particle_charge = np.concatenate((_particle_charge, _data[i].array('Particle.Charge')))

                    _jet_event = np.concatenate((_jet_event,_data[i].array('Event')))
                    _jet_pt = np.concatenate((_jet_pt, _data[i].array('Jet.PT')))
                    _jet_eta = np.concatenate((_jet_eta, _data[i].array('Jet.Eta')))
                    _jet_phi = np.concatenate((_jet_phi, _data[i].array('Jet.Phi')))
                    _jet_btag = np.concatenate((_jet_btag, _data[i].array('Jet.BTag')))
                    _jet_area = np.concatenate((_jet_area, _data[i].array('Jet.Area')))
                    _jet_mass = np.concatenate((_jet_mass, _data[i].array('Jet.Mass')))
                    _jet_charge = np.concatenate((_jet_charge, _data[i].array('Jet.Charge')))
                    _num_of_jets =  np.concatenate((_num_of_jets, _data[i].array('Jet')))
                
            if self.model in self.require_lepton:
                for i in tqdm.trange(num_of_files):
                    if i == 0 :
                        _missing_et_met = _data[i].array('MissingET.MET')
                        _missing_et_eta = _data[i].array('MissingET.Eta')
                        _missing_et_phi = _data[i].array('MissingET.Phi')
                        
                        _muon_pt = _data[i].array('Muon.PT')
                        _muon_eta = _data[i].array('Muon.Eta')
                        _muon_phi = _data[i].array('Muon.Phi')
                        _muon_charge = _data[i].array('Muon.Charge')        
    
                        _electron_pt = _data[i].array('Electron.PT')
                        _electron_eta = _data[i].array('Electron.Eta')
                        _electron_phi = _data[i].array('Electron.Phi')
                        _electron_charge = _data[i].array('Electron.Charge')
                    else: 
                        _missing_et_met = np.concatenate((_missing_et_met, _data[i].array('MissingET.MET')))
                        _missing_et_eta = np.concatenate((_missing_et_eta, _data[i].array('MissingET.Eta')))
                        _missing_et_phi = np.concatenate((_missing_et_phi, _data[i].array('MissingET.Phi')))
                        
                        _muon_pt = np.concatenate((_muon_pt, _data[i].array('Muon.PT')))
                        _muon_eta = np.concatenate((_muon_eta, _data[i].array('Muon.Eta')))
                        _muon_phi = np.concatenate((_muon_phi, _data[i].array('Muon.Phi')))
                        _muon_charge = np.concatenate((_muon_charge, _data[i].array('Muon.Charge')))

                        _electron_pt = np.concatenate((_electron_pt, _data[i].array('Electron.PT')))
                        _electron_eta = np.concatenate((_electron_eta, _data[i].array('Electron.Eta')))
                        _electron_phi = np.concatenate((_electron_phi, _data[i].array('Electron.Phi')))
                        _electron_charge = np.concatenate((_electron_charge, _data[i].array('Electron.Charge')))
        else:
            _data = uproot.open(self.path)['Delphes']
            
            print("Loading particle information.")
            _particle_event = _data.array('Event')
            _particle_pt = _data.array('Particle.PT')
            _particle_eta = _data.array('Particle.Eta')
            _particle_phi = _data.array('Particle.Phi')
            _particle_pid = _data.array('Particle.PID')
            _particle_M1 = _data.array('Particle.M1')
            _particle_M2 = _data.array('Particle.M2')
            _particle_D1 = _data.array('Particle.D1')
            _particle_D2 = _data.array('Particle.D2')
            _particle_status = _data.array('Particle.Status')
            _particle_rapidity = _data.array('Particle.Rapidity')
            _particle_mass = _data.array('Particle.Mass')
            _particle_charge = _data.array('Particle.Charge')
            
            print("Loading jet information.")
            _jet_event = _data.array('Event')
            _jet_pt = _data.array('Jet.PT')
            _jet_eta = _data.array('Jet.Eta')
            _jet_phi = _data.array('Jet.Phi')
            _jet_btag = _data.array('Jet.BTag')
            _jet_area = _data.array('Jet.Area')
            _jet_mass = _data.array('Jet.Mass')
            _jet_charge = _data.array('Jet.Charge')
            _num_of_jets = _data.array('Jet')
            
            if self.model in self.require_lepton:

                print("Loading MET information.")
                _missing_et_met = _data.array('MissingET.MET')
                _missing_et_eta = _data.array('MissingET.Eta')
                _missing_et_phi = _data.array('MissingET.Phi')

                print("Loading muon information.")
                _muon_pt = _data.array('Muon.PT')
                _muon_eta = _data.array('Muon.Eta')
                _muon_phi = _data.array('Muon.Phi')
                _muon_charge = _data.array('Muon.Charge')

                print("Loading electron information.")
                _electron_pt = _data.array('Electron.PT')
                _electron_eta = _data.array('Electron.Eta')
                _electron_phi = _data.array('Electron.Phi')
                _electron_charge = _data.array('Electron.Charge')        

        jet_dataset = {
            "event": _jet_event, 
            "pt": _jet_pt, 
            "eta": _jet_eta, 
            "phi": _jet_phi, 
            "btag": _jet_btag, 
            "area": _jet_area, 
            "mass": _jet_mass, 
            "num_of_jets": _num_of_jets, 
            "charge": _jet_charge,                                     
        }
        particle_dataset = {
            "event": _particle_event, 
            "pt": _particle_pt, 
            "eta": _particle_eta, 
            "phi": _particle_phi, 
            "pid": _particle_pid, 
            "M1": _particle_M1, 
            "M2": _particle_M2, 
            "D1": _particle_D1, 
            "D2": _particle_D2, 
            "status": _particle_status, 
            "rapidity": _particle_rapidity, 
            "mass": _particle_mass, 
            "charge": _particle_charge,                                      
        }
        if self.model in self.require_lepton:
            muon_dataset = {
                "pt": _muon_pt, 
                "eta": _muon_eta, 
                "phi": _muon_phi, 
                "charge": _muon_charge,
            }
            electron_dataset = {
                "pt": _electron_pt, 
                "eta": _electron_eta, 
                "phi": _electron_phi, 
                "charge": _electron_charge,
            }
            MissingET_dataset = {
                "MET": _missing_et_met, 
                "eta": _missing_et_eta, 
                "phi": _missing_et_phi, 
            }
            dataset = {
                "particle": particle_dataset,
                "jet": jet_dataset,
                "muon": muon_dataset,
                "electron": electron_dataset,
                "MissingET": MissingET_dataset,
            }
        else: 
            dataset = {
                "particle": particle_dataset,
                "jet": jet_dataset,
            }
        return dataset
    def wrtite_hdf5(self):
        if self.model == 'ttbar':
            jet_index_dict_name = OrderedDict((
                ("left_target", ["b", "q1", "q2"]),
                ("right_target", ["b", "q1", "q2"]),
            ))
        elif self.model == 'ttH':
            jet_index_dict_name = OrderedDict((
                ("left_target", ["b", "q1", "q2"]),
                ("right_target", ["b", "q1", "q2"]),
                ("higgs_target", ["b1", "b2"]),
            ))
        elif self.model == 'four_top':
            jet_index_dict_name = OrderedDict((
                ("first_target", ["b", "q1", "q2"]),
                ("second_target", ["b", "q1", "q2"]),
                ("third_target", ["b", "q1", "q2"]),
                ("fourth_target", ["b", "q1", "q2"]),
            ))
        elif self.model == 'ttbar_lep':
            jet_index_dict_name = OrderedDict((
                ("left_target", ["b", "neutrino", "lepton"]),
                ("right_target", ["b", "q1", "q2"]),
            ))
        elif self.model == 'ttbar_lep_left':
            jet_index_dict_name = OrderedDict((
                ("left_target", ["b", "neutrino", "lepton"]),
                ("right_target", ["b", "q1", "q2"]),
            ))
        elif self.model == 'ttbar_lep_right':
            jet_index_dict_name = OrderedDict((
                ("left_target", ["b", "q1", "q2"]),
                ("right_target", ["b", "lepton", "neutrino"]),
            ))
        else:
            print("Please select a correct model.")
        parton_features = self.kargs['parton_features']
        jet_features = self.kargs['jet_features']
        if self.model in self.require_lepton:
            lepton_features = self.kargs['lepton_features']
            met_features = self.kargs['MET_features']
        if self.kargs['usage'] == 'parse':
            targets = self.kargs['targets']
            masks = self.kargs['masks']
            if self.kargs['contain_chi_square']:
                chi_square_result = self.kargs['chi_square_result']
        elif self.kargs['usage'] == 'chi2':
            if self.kargs['contain_chi_square'] == False:
                pass
            else:
                chi_square_result = self.kargs['chi_square_result']
        else: 
            print("Please input a correct parameter for writing record to hdf5 files.")
        with h5py.File(self.path, 'w') as file:
            for key, value in parton_features.items():
                file.create_dataset(f"parton_features/{key}", data=value)
            for key, value in jet_features.items():
                file.create_dataset(f"jet_features/{key}", data=value)  
            if self.model in self.require_lepton:
                for key, value in lepton_features.items():
                    file.create_dataset(f"lepton_features/{key}", data=value)  
                for key, value in met_features.items():
                    file.create_dataset(f"met_features/{key}", data=value)  
            if self.kargs['usage'] == 'parse':
                for key, value in jet_index_dict_name.items():
                    target = targets[key]
                    mask = masks[key]
                    file.create_dataset(f"target/{key}/mask", data=mask)
                    for quark_name, jet in zip(value, target):
                        file.create_dataset(f"target/{key}/{quark_name}", data=jet)
                if self.kargs['contain_chi_square']:
                    for key, value in chi_square_result.items():
                        file.create_dataset(f"chi_square_result/{key}", data=value) 
            if self.kargs['usage'] == 'chi2' and self.kargs['contain_chi_square']:
                for key, value in chi_square_result.items():
                    file.create_dataset(f"chi_square_result/{key}", data=value) 

@nb.jit(nopython=True)
def deltaPhi(phi1: float,phi2: float) -> float:
    """
    This is a function deltaPhi value between two target.
    phi1: phi value from target one. 
    phi2: phi value from target two. 
    """
    phi = phi1-phi2
    while phi >= np.pi: phi -= np.pi*2.
    while phi < -np.pi: phi += np.pi*2.
    return phi

@nb.jit(nopython=True)
def delta_R(eta1: float, phi1: float, eta2: float, phi2: float) -> float:
    """
    This is a function delta_R value between two target.
    phi1: phi value from target one. 
    eta1: eta value from target one. 
    phi2: phi value from target two. 
    eta2: eta value from target two. 
    """
    return float(np.sqrt(deltaPhi(phi1,phi2)**2+(eta1-eta2)**2))

class process_methods:
    @nb.jit(nopython=True)
    def __jet_marker(Pt: np.ndarray, Eta: np.ndarray, Btag: np.ndarray) -> Union[np.ndarray, np.ndarray]:
        _marker_jet = []
        _marker_btag = []

        for j in range(len(Pt)):
            if Btag[j] == 1 and Pt[j] > 25 and np.abs(Eta[j]) < 2.5:
                _marker_btag.append(1) 
            else: 
                _marker_btag.append(0) 

            if Pt[j] > 25 and np.abs(Eta[j]) <= 2.5:
                _marker_jet.append(1)
            else:
                _marker_jet.append(0)
        return np.array(_marker_jet), np.array(_marker_btag)
    @nb.jit(nopython=True)
    def __event_marker(_jet: np.ndarray, _btag: np.ndarray, _requirement1: int, _requirement2: int) -> int:
        if np.sum(_jet == 1) >= _requirement1 and np.sum(_btag == 1) >= _requirement2 :
            _result = 1
        else:
            _result = 0
        return _result
    
    @staticmethod
    def event_selection(MODEL: str, **kargs) -> np.ndarray:
        marker_event = []
        marker_jet = []
        marker_btag = []
        PT = kargs['pt']
        ETA = kargs['eta']
        BTAG = kargs['btag']
        print("MODE: {0}, Number of events: {1}.".format(MODEL, len(PT)))
        
        requirement = {
            "ttbar": [2, 6],
            "ttH": [2, 8],
            "four_top": [2, 12],
            "ttbar_lep": [2, 4],
            "ttbar_lep_left": [2, 4],
            "ttbar_lep_right": [2, 4], 
            "ZH": [2, 6],
        }     
        
        if MODEL != 'ttbar_lep' and MODEL != 'ttbar_lep_left' and MODEL != "ttbar_lep_right":
            for i in tqdm.trange(len(PT), desc="Marking jets"):
                try:
                    __marker_jet, __marker_bjet = process_methods.__jet_marker(np.array(PT[i]), np.array(ETA[i]), np.array(BTAG[i]))
                except:
                    # There is some empty list exist, use this exception to avoid error.
                    __marker_jet = np.zeros(6)
                    __marker_jet = np.zeros(6)
                    
                marker_jet.append(__marker_jet)
                marker_btag.append(__marker_bjet)

            marker_jet = np.asanyarray(marker_jet, dtype='object')
            marker_btag = np.asanyarray(marker_btag, dtype='object')
            
            for i in tqdm.trange(len(PT), desc="Marking events"):
                marker_event.append(process_methods.__event_marker(marker_jet[i], marker_btag[i], requirement[MODEL][1], requirement[MODEL][0]))
            marker_event = np.asanyarray(marker_event, dtype=object)
        else: 
            PHI = kargs['phi']

            ELECTRON_PT = [ list(x) if x.size>=2 else x.item() if x.size==1 else 99999 for x in kargs['electron_pt']]
            ELECTRON_ETA = [ list(x) if x.size>=2 else x.item() if x.size==1 else 99999 for x in kargs['electron_eta']]
            ELECTRON_PHI = [ list(x) if x.size>=2 else x.item() if x.size==1 else 99999 for x in kargs['electron_phi']]

            MUON_PT = [ list(x) if x.size>=2 else x.item() if x.size==1 else 99999 for x in kargs['muon_pt']]
            MUON_ETA = [ list(x) if x.size>=2 else x.item() if x.size==1 else 99999 for x in kargs['muon_eta']]
            MUON_PHI = [ list(x) if x.size>=2 else x.item() if x.size==1 else 99999 for x in kargs['muon_phi']]

            marker_lepton = []
            LEPTON_PT = []
            LEPTON_ETA = []
            LEPTON_PHI = []
            for a,b in zip(ELECTRON_PT, MUON_PT):
                _tmp = []
                if isinstance(a, float) or isinstance(a, int):
                    _tmp.append(a)
                elif isinstance(a, list):
                    for c in a:
                        _tmp.append(c)
                else: 
                    print('error', type(a))

                if isinstance(b, float) or isinstance(b, int):
                    _tmp.append(b)
                elif isinstance(b, list):
                    for c in b:
                        _tmp.append(c)
                else: 
                    print('error', type(b))
                LEPTON_PT.append(_tmp)

            for a,b in zip(ELECTRON_ETA, MUON_ETA):
                _tmp = []
                if isinstance(a, float) or isinstance(a, int):
                    _tmp.append(a)
                elif isinstance(a, list):
                    for c in a:
                        _tmp.append(c)
                else: 
                    print('error', type(a))

                if isinstance(b, float) or isinstance(b, int):
                    _tmp.append(b)
                elif isinstance(b, list):
                    for c in b:
                        _tmp.append(c)
                else: 
                    print('error', type(b))
                LEPTON_ETA.append(_tmp)
            for a,b in zip(ELECTRON_PHI, MUON_PHI):
                _tmp = []
                if isinstance(a, float) or isinstance(a, int):
                    _tmp.append(a)
                elif isinstance(a, list):
                    for c in a:
                        _tmp.append(c)
                else: 
                    print('error', type(a))

                if isinstance(b, float) or isinstance(b, int):
                    _tmp.append(b)
                elif isinstance(b, list):
                    for c in b:
                        _tmp.append(c)
                else: 
                    print('error', type(b))
                LEPTON_PHI.append(_tmp)
            for i in tqdm.trange(len(PT), desc="Marking jets"):
                try:
                    __marker_jet, __marker_bjet = process_methods.__jet_marker(np.array(PT[i]), np.array(ETA[i]), np.array(BTAG[i]))
                    
                except:
                    __marker_jet = np.zeros(6)
                    __marker_jet = np.zeros(6)
                marker_jet.append(__marker_jet)
                marker_btag.append(__marker_bjet)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_btag = np.asanyarray(marker_btag, dtype=object)

            #Remove electron from jets catogary
            for i in tqdm.trange(len(PT), desc='Doing lepton-jet removal'):
                for j in range(len(PT[i])):
                    for k in range(len(LEPTON_PT[i])):
                        if delta_R(ETA[i][j], PHI[i][j], LEPTON_ETA[i][k], LEPTON_PHI[i][k]) < 0.4:
                            marker_jet[i][j] = 0
                        else : pass 

            for i in tqdm.trange(len(LEPTON_PT), desc='Marking leptons'):
                _marker_lepton = []
                for j in range(len(LEPTON_PT[i])):
                    if LEPTON_PT[i][j] > 25 and np.abs(LEPTON_ETA[i][j]) < 2.5:
                        _marker_lepton.append(1)
                    else :
                        _marker_lepton.append(0)
                marker_lepton.append(np.asanyarray(_marker_lepton, dtype=object))
            marker_lepton = np.asanyarray(marker_lepton, dtype=object)
            for i in tqdm.trange(len(PT), desc="Marking events"):
                if np.sum(marker_jet[i] == 1) >= requirement[MODEL][1] and np.sum(marker_btag[i] == 1) >= requirement[MODEL][0] and np.sum(marker_lepton[i] ==1) == 1 and (np.array(LEPTON_PT[i]) != 99999).sum() == 1:
                    marker_event.append(1)
                else:
                    marker_event.append(0)
            marker_event = np.asanyarray(marker_event, dtype=object)
        return marker_event
    
    @staticmethod
    def shifted_particle_tracing(dataset: pd.core.frame.DataFrame, PID_daughter: int, idx: int) -> int:
        """
        This is a frunction tracing the on-flying particle. 
        """
        if (dataset.iloc[idx,6] == PID_daughter):
            return int(dataset.iloc[idx,4])
        
    @staticmethod
    def daughter_finder(dataset: pd.core.frame.DataFrame, PID: int, MODEL: str, not_stable_FSP=True) -> dict:
        if not_stable_FSP:
            if MODEL == 'four_top':
                mother_cand_idx = np.array(dataset[dataset["PID"] == PID]["Index"])[-10:]
                mother_idx = []
                filtered_idx = []
                for _idx, element in enumerate(mother_cand_idx):
                    idx_1 = dataset['Daughter_1'][element]
                    idx_2 = dataset['Daughter_2'][element]
                    pid_1 = dataset["PID"].loc[idx_1]
                    pid_2 = dataset["PID"].loc[idx_2]
                    blacklist = [21, 22]
                    if idx_1 != idx_2 and pid_1 not in blacklist and pid_2 not in blacklist:
                        filtered_idx.append(_idx)
                        mother_idx.append(int(element))
                D1_1 = dataset["Daughter_1"].iloc[mother_idx[0]]
                D1_2 = dataset["Daughter_2"].iloc[mother_idx[0]]
                D2_1 = dataset["Daughter_1"].iloc[mother_idx[1]]
                D2_2 = dataset["Daughter_2"].iloc[mother_idx[1]]
                if dataset["PID"][int(D1_1)] < dataset["PID"][int(D1_2)]:
                    D1_1, D1_2 = D1_2, D1_1
                if dataset["PID"][int(D2_1)] < dataset["PID"][int(D2_2)]:
                    D2_1, D2_2 = D2_2, D2_1
                _result = {
                    "mother_1_idx": mother_idx[0],
                    "daughter_1_1_idx": D1_1,
                    "daughter_1_2_idx": D1_2,
                    "mother_2_idx": mother_idx[1],
                    "daughter_2_1_idx": D2_1,
                    "daughter_2_2_idx": D2_2,
                }
            else:
                mother_idx = np.array(dataset[dataset["PID"] == PID]["Index"])[-1]
                D1 = np.array(dataset[dataset["PID"] == PID]["Daughter_1"])[-1]
                D2 = np.array(dataset[dataset["PID"] == PID]["Daughter_2"])[-1]
                if dataset["PID"][int(D1)] < dataset["PID"][int(D2)]:
                    D1, D2 = D2, D1
                _result = {
                    "mother_idx": mother_idx,
                    "daughter_1_idx": D1,
                    "daughter_2_idx": D2,
                }
        else:
            STATUS = kargs['status']
            MODEL = kargs['model']
            if MODEL == "ttbar_lep_right" or MODEL == "ttbar_lep_left":
                for i in range(len(dataset)):
                    if(dataset.iloc[i,1] == STATUS and dataset.iloc[i,6] == PID ): 
                        daughter_index = int(dataset.iloc[i,0])
                if( dataset.iloc[daughter_index,6] == PID ):
                    shifted_particle_index = dataset.iloc[daughter_index, 4]


                while dataset.iloc[shifted_particle_index,6] == PID:
                    init_shifted_particle_index = shifted_particle_index
                    shifted_particle_index = shifted_particle_tracing(dataset, PID, init_shifted_particle_index)       

                dauthter_idx_1 = dataset.iloc[init_shifted_particle_index, 4]
                daughter_pid_1 = dataset.iloc[dauthter_idx_1, 6]

                dauthter_idx_2 = dataset.iloc[init_shifted_particle_index, 5]
                daughter_pid_2 = dataset.iloc[dauthter_idx_2, 6]

                _result = {
                    "mother_idx": init_shifted_particle_index,
                    "daughter_1_idx": dauthter_idx_1,
                    "daughter_2_idx": dauthter_idx_2,
                }
        return _result

    @staticmethod
    def lephad_finder(dataset: pd.core.frame.DataFrame, PID: int, LEPHAD: int, MODEL: str, not_stable_FSP=True) -> dict:

    ## PID: only accepts top and W here. The sign does not matter. 
    ## LEPHAD: 1 for leptonic, 2 for hadronic 

        # Extract W+ info 
        mother_idx = np.array(dataset[dataset["PID"] == 24]["Index"])[-1]
        D1 = np.array(dataset[dataset["PID"] == 24]["Daughter_1"])[-1]
        D2 = np.array(dataset[dataset["PID"] == 24]["Daughter_2"])[-1]
        # Extract W- info
        mother_idx_m = np.array(dataset[dataset["PID"] == -24]["Index"])[-1]
        D1m = np.array(dataset[dataset["PID"] == -24]["Daughter_1"])[-1]
        D2m = np.array(dataset[dataset["PID"] == -24]["Daughter_2"])[-1]
        
        #print("W+ daughter 1 PID", dataset["PID"][int(D1)])
        #print("W+ daughter 2 PID", dataset["PID"][int(D2)])
        #print("W- daughter 1 PID", dataset["PID"][int(D1m)])
        #print("W- daughter 2 PID", dataset["PID"][int(D2m)])

	# Extract top info 
        mother_idx_t = np.array(dataset[dataset["PID"] == 6]["Index"])[-1]
        D1t = np.array(dataset[dataset["PID"] == 6]["Daughter_1"])[-1]
        D2t = np.array(dataset[dataset["PID"] == 6]["Daughter_2"])[-1]
        # Extract anti-top info
        mother_idx_at = np.array(dataset[dataset["PID"] == -6]["Index"])[-1]
        D1at = np.array(dataset[dataset["PID"] == -6]["Daughter_1"])[-1]
        D2at = np.array(dataset[dataset["PID"] == -6]["Daughter_2"])[-1]

        #print("t daughter 1 PID", dataset["PID"][int(D1t)])
        #print("t daughter 2 PID", dataset["PID"][int(D2t)])
        #print("t~ daughter 1 PID", dataset["PID"][int(D1at)])
        #print("t~ daughter 2 PID", dataset["PID"][int(D2at)])

        if abs(PID)==24: 
            if LEPHAD ==1 and abs(dataset["PID"][int(D1)]) > 10 : #record leptonic W & its daughters
                if abs(dataset["PID"][int(D1)]) < abs(dataset["PID"][int(D2)]):
                    D1, D2 = D2, D1
                _result = {
                    "mother_idx": mother_idx,
                    "daughter_1_idx": D1,
                    "daughter_2_idx": D2,
                }

            elif LEPHAD ==1 and abs(dataset["PID"][int(D1m)]) > 10 : #record leptonic W & its daughters
                if abs(dataset["PID"][int(D1m)]) < abs(dataset["PID"][int(D2m)]):
                    D1m, D2m = D2m, D1m
                _result = {
                    "mother_idx": mother_idx_m,
                    "daughter_1_idx": D1m,
                    "daughter_2_idx": D2m,
                }

            elif LEPHAD ==2 and abs(dataset["PID"][int(D1)]) < 10 : #record hadronic W & its daughters
                if abs(dataset["PID"][int(D1)]) < abs(dataset["PID"][int(D2)]):
                    D1, D2 = D2, D1
                _result = {
                    "mother_idx": mother_idx,
                    "daughter_1_idx": D1,
                    "daughter_2_idx": D2,
                }

            elif LEPHAD ==2 and abs(dataset["PID"][int(D1m)]) < 10 : #record hadronic W & its daughters
                if abs(dataset["PID"][int(D1m)]) < abs(dataset["PID"][int(D2m)]):
                    D1m, D2m = D2m, D1m
                _result = {
                    "mother_idx": mother_idx_m,
                    "daughter_1_idx": D1m,
                    "daughter_2_idx": D2m,
                }

        elif abs(PID)==6: 
            if LEPHAD ==1 and abs(dataset["PID"][int(D1)]) > 10 : #record leptonic top, using D1 to check W+ decay
                if abs(dataset["PID"][int(D1t)]) < abs(dataset["PID"][int(D2t)]):
                    D1t, D2t = D2t, D1t
                _result = {
                    "mother_idx": mother_idx_t,
                    "daughter_1_idx": D1t,
                    "daughter_2_idx": D2t,
                }   

            elif LEPHAD ==1 and abs(dataset["PID"][int(D1m)]) > 10 : #record leptonic top, using D1m to check W- decay
                if abs(dataset["PID"][int(D1at)]) < abs(dataset["PID"][int(D2at)]):
                    D1at, D2at = D2at, D1at
                _result = {
                    "mother_idx": mother_idx_at,
                    "daughter_1_idx": D1at,
                    "daughter_2_idx": D2at,
                }

            elif LEPHAD ==2 and abs(dataset["PID"][int(D1)]) < 10 : #record hadronic top, using D1 to check W+ decay
                if abs(dataset["PID"][int(D1t)]) < abs(dataset["PID"][int(D2t)]):
                    D1t, D2t = D2t, D1t
                _result = {
                    "mother_idx": mother_idx_t,
                    "daughter_1_idx": D1t,
                    "daughter_2_idx": D2t,
                }

            elif LEPHAD ==2 and abs(dataset["PID"][int(D1m)]) < 10 : #record hadronic top, using D1m to check W- decay
                if abs(dataset["PID"][int(D1at)]) < abs(dataset["PID"][int(D2at)]):
                    D1at, D2at = D2at, D1at
                _result = {
                    "mother_idx": mother_idx_at,
                    "daughter_1_idx": D1at,
                    "daughter_2_idx": D2at,
                }

        return _result


    #tracing the daughters
    #Input two daughter of top/top_bar and find their daughter
    @staticmethod
    def quark_finder(dataset: pd.core.frame.DataFrame, mother_idx_1: int, mother_idx_2: int) -> dict:
        """
        This is a function finding the daughter of bosons.
        This function will no longer be deprecated. 
        Check before calling this function.
        """
        #Specific two daughter of top
        def W_b_specifier(dataset, input_1_idx, input_2_idx):
            if dataset.iloc[int(input_1_idx),6] == PID.w_plus or dataset.iloc[int(input_1_idx),6] == PID.w_minus :
                return int(input_1_idx), int(dataset.iloc[int(input_1_idx),6]), int(input_2_idx)
            elif dataset.iloc[int(input_1_idx),6] == PID.bottom or dataset.iloc[int(input_1_idx),6] == PID.anti_bottom :
                return  int(input_2_idx), int(dataset.iloc[int(input_1_idx),6]), int(input_1_idx)
            else :
                pass
                #print("Please check your data.")

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
        _result = {
                "b_idx": b_quark_idx,
                "W_daughter_1_idx": daughter_1_idx,
                "W_dauthter_2_idx": daughter_2_idx,
            }
        return  _result

    
class helper:
    @staticmethod
    def to_dataframe(DATASET: dict, index: int) -> pd.core.frame.DataFrame:
        idx = np.linspace(0, len( DATASET["pt"][index])-1, num = len( DATASET["pt"][index]) )
        dic = {
             "Index": idx,
            "Status":  DATASET["status"][index],
            "Mother_1":  DATASET["M1"][index],
            "Mother_2":  DATASET["M2"][index],
            "Daughter_1":  DATASET["D1"][index],
            "Daughter_2":  DATASET["D2"][index],
            "PID": DATASET["pid"][index],
            "PT":  DATASET["pt"][index],
            "Eta":  DATASET["eta"][index],
            "Phi":  DATASET["phi"][index],
            "Mass":  DATASET["mass"][index],
        }
        return pd.DataFrame(dic)
    @staticmethod
    def fetch_kinematics_properties_from_dataset(dataset: pd.core.frame.DataFrame, index: np.ndarray, colume_name: str)-> list:
        _result = [dataset[colume_name][_index] for _index in index]
        return _result

def old_deltaR_matching(NUM_OF_PARTON, NUM_OF_JET, PARTON_ETA, PARTON_PHI, JET_ETA, JET_PHI, CUTS, MODEL):
    """
    This is a function for doing delta R matching.
    PARTON_ETA: Array, a list of partons's eta in a event.
    PARTON_PHI: Array, a list of partons's phi in a event.
    JET_ETA: Array, a list of jet's eta in a event.
    JET_PHI: Array, a list of jet's phi in a event.
    """
    _dR_between_parton_jet = []
    
    _parton_jet_index = np.full(NUM_OF_PARTON, -1)
    _jet_parton_index = np.full(NUM_OF_JET, -1)
    
    _jet_to_parton_list = np.zeros(len(PARTON_ETA))
    _parton_to_jet_list = np.zeros(len(JET_ETA))
    _min_value = []
    j = 0
    a = 0
    b = 0
    while a < NUM_OF_PARTON :
        for b in range( NUM_OF_JET ):
            _dR_between_parton_jet.append(calc_helper.delta_R( PARTON_ETA[a], PARTON_PHI[a], JET_ETA[b], JET_PHI[b]))
            j +=1
        a += 1 
    
    array = np.reshape(np.array(_dR_between_parton_jet), [NUM_OF_PARTON, NUM_OF_JET])
    array_index = [x for x in range(len(PARTON_ETA))]

    _dataset = pd.DataFrame(index = array_index, data = array).T
    for j in range(len(PARTON_ETA)):
        min_val = _dataset.stack().min()
        if min_val < CUTS:
            min_idx, min_col = _dataset.stack().idxmin()
            
            _jet_to_parton_list[j] = int(min_idx)
            _parton_to_jet_list[j] = int(min_col)
            _dataset = _dataset.drop([min_col], axis=1)
            _dataset = _dataset.drop([min_idx], axis=0)
            _min_value.append(min_val)
        else:
            _jet_to_parton_list[j] = 'Nan'
            _parton_to_jet_list[j] = 'Nan'
    for k in range(NUM_OF_PARTON, NUM_OF_JET):
        _parton_to_jet_list[k] = 'Nan'
    
    if MODEL != 'ttbar_lep' and MODEL != 'ttbar_lep_left' and MODEL != 'ttbar_lep_right':
        for j in range(len(JET_ETA)):
            for k in range(NUM_OF_JET):
                if _parton_to_jet_list[j] == k:
                    _parton_jet_index[k] = int(_jet_to_parton_list[j])                    
    elif MODEL == 'ttbar_lep' or MODEL == 'ttbar_lep_left' or MODEL == 'ttbar_lep_right':
        for j in range(len(JET_ETA)):
            if _parton_to_jet_list[j] == 0 :
                _parton_jet_index[0] = int(_jet_to_parton_list[j])
            else: 
                pass

            if _parton_to_jet_list[j] == 1 :
                _parton_jet_index[1] = int(_jet_to_parton_list[j])
            else: 
                pass
            if _parton_to_jet_list[j] == 2 :
                _parton_jet_index[2] = int(_jet_to_parton_list[j])
            else: 
                pass

            if _parton_to_jet_list[j] == 3 :
                _parton_jet_index[3] = int(_jet_to_parton_list[j])
            else:
                pass
    else:
        print("Delta R matching faild, please check your model.")

    ll = len(JET_ETA)
    for k in range(NUM_OF_PARTON):
        for m in range(ll):
            if _jet_to_parton_list[k] == int(m):
                _jet_parton_index[int(m)] = _parton_to_jet_list[k]
            else: pass

    
    return _jet_parton_index, _parton_jet_index

@nb.jit(nopython=True)
def deltaR_matching(NUM_OF_DAUGHTER, NUM_OF_JET, PARTON_ETA, PARTON_PHI, PARTON_PID, JET_ETA, JET_PHI, CUTS, MODEL):
    """
    This is a function for doing delta R matching.
    PARTON_ETA: Array, a list of partons's eta in a event.
    PARTON_PHI: Array, a list of partons's phi in a event.
    JET_ETA: Array, a list of jet's eta in a event.
    JET_PHI: Array, a list of jet's phi in a event.
    """
    _dR_between_parton_jet = []   
    _parton_jet_index = np.full(NUM_OF_DAUGHTER, -1)
    _jet_parton_index = np.full(NUM_OF_JET, -1)

    _jet_to_parton_list = np.zeros(len(PARTON_ETA))
    _parton_to_jet_list = np.zeros(len(JET_ETA))
    _min_value = []
    j = 0
    a = 0
    b = 0
    while a < len(PARTON_ETA) :
        for b in range( NUM_OF_JET ):
            if abs(PARTON_PID[a]) < 6 :
               _dR_between_parton_jet.append(delta_R(PARTON_ETA[a], PARTON_PHI[a], JET_ETA[b], JET_PHI[b]))
            else :
               _dR_between_parton_jet.append(10.) # exclude leptons by assigning large dR
            j +=1
        a += 1 
    src_array = np.array(_dR_between_parton_jet)
    min_pos_list = np.argsort(src_array)
    _parton_to_jet_list = []
    _jet_to_parton_list = []
    for i in range(len(min_pos_list)):
        _x = min_pos_list[i]//NUM_OF_JET
        _y = min_pos_list[i]%NUM_OF_JET
        if _x not in _jet_to_parton_list and _y not in _parton_to_jet_list and src_array[min_pos_list[i]] < CUTS:
            _jet_to_parton_list.append(_x)
            _parton_to_jet_list.append(_y)
            _min_value.append(src_array[min_pos_list[i]])
    for a, b in zip(_parton_to_jet_list, _jet_to_parton_list):
        pair_jet = a
        pair_parton = b
        _parton_jet_index[b] = a
        _jet_parton_index[a] = b

    return _jet_parton_index, _parton_jet_index

def gen_combinaton(NUM_OF_JETS, MODEL, BTAG_JET):
    num_of_btag = np.sum(np.array(BTAG_JET) ==1)
    if MODEL == 'ttbar':
        _bjet_list = []  
        _jet_list = []
        bjet = []
        _jet_index = []
        
        for i in range(NUM_OF_JETS):
            _jet_index.append(i)

        for i in range(NUM_OF_JETS):
            if BTAG_JET[i] == 1:
                _bjet_list.append(i)
            else :
                _jet_list.append(i)
                
        _bjet = itertools.combinations(_bjet_list, 2)
        
        for a in _bjet:
            bjet.append(a)

        bjet = np.array(bjet, dtype='object')

        _jet_index_candidate = []

        for i in range(len(bjet)):

            jet = []
            
            tmp_jet_index = _bjet_list.copy()
            for c in range(len(bjet[i])): 
                
                _tmp = bjet[i][c]
                tmp_jet_index.remove(_tmp)
            
            tmp_jet_list = _jet_list.copy()

            for d in tmp_jet_index:
                tmp_jet_list.append(d)
            
            _jet = itertools.permutations(tmp_jet_list, 4)

            for b in _jet:
                jet.append(b)

            jet = np.array(jet, dtype='object')

            for j in range(len(jet)):
                tmp_jet_index_candidate = []
                tmp_jet_index_candidate.append(bjet[i][0])
                tmp_jet_index_candidate.append(jet[j][0])
                tmp_jet_index_candidate.append(jet[j][1])
                tmp_jet_index_candidate.append(bjet[i][1])
                tmp_jet_index_candidate.append(jet[j][2])
                tmp_jet_index_candidate.append(jet[j][3])
                _jet_index_candidate.append(np.array(tmp_jet_index_candidate))
        _jet_index_candidate = np.array(_jet_index_candidate)
    elif MODEL == 'ttH':
        tmp_jet_list = []
        tmp_bjet_list = []
        
        _jet_index = []

        for i in range(NUM_OF_JETS):
            _jet_index.append(i)

        for i in range(NUM_OF_JETS):
            if BTAG_JET[i] == 1:
                tmp_bjet_list.append(i)
            else :
                tmp_jet_list.append(i)
        
        if num_of_btag == 4:
            b_jet = []
            jet = []
            _bjet = itertools.combinations(tmp_bjet_list, 4)
            _jet = itertools.combinations(tmp_jet_list, 4)
            
            b_jet = np.array([x for x in _bjet], dtype='object')
            jet = np.array([x for x in _jet], dtype='object')

            _jet_index_candidate = []
            for i in range(len(b_jet)):
                for j in range(len((jet))):
                    tmp_jet_index_candidate = []
                    tmp_jet_index_candidate.append(b_jet[i][0])
                    tmp_jet_index_candidate.append(jet[j][0])
                    tmp_jet_index_candidate.append(jet[j][1])
                    tmp_jet_index_candidate.append(b_jet[i][1])
                    tmp_jet_index_candidate.append(jet[j][2])
                    tmp_jet_index_candidate.append(jet[j][3])
                    tmp_jet_index_candidate.append(b_jet[i][2])
                    tmp_jet_index_candidate.append(b_jet[i][3])
                    _jet_index_candidate.append(np.array(tmp_jet_index_candidate))
            _jet_index_candidate = np.array(_jet_index_candidate)
        elif num_of_btag != 4:
            
            require_num_bjet = 4
            lack_of_bjet = require_num_bjet - num_of_btag 
            if lack_of_bjet > 0:
                _jet_for_append = itertools.combinations(tmp_jet_list, lack_of_bjet)

                jet_for_append = np.array([x for x in _jet_for_append], dtype='object')
                for i in range(len(jet_for_append)):
                    b_jet = []
                    jet = []
                    tmp_jet_index = tmp_jet_list.copy()
                    for c in range(len(jet_for_append[i])): 
                        _tmp = jet_for_append[i][c]
                        tmp_jet_index.remove(_tmp)
                    
                    tmp_bjet_index = tmp_bjet_list.copy()

                    for d in jet_for_append[i]:                    
                        tmp_bjet_index.append(int(d))
                    

                    _bjet = itertools.permutations(tmp_bjet_index, 4)
                    _jet = itertools.permutations(tmp_jet_index, 4)
                    
                    for a, b in zip(_jet, _bjet):
                        jet.append(a)
                        b_jet.append(b)
                    
                    jet = np.array(jet, dtype='object')
                    b_jet = np.array(b_jet, dtype='object')
                    
                    _jet_index_candidate = []
                    for j in range(len(b_jet)):
                        for k in range(len((jet))):
                            tmp_jet_index_candidate = []
                            tmp_jet_index_candidate.append(b_jet[j][0])
                            tmp_jet_index_candidate.append(jet[k][0])
                            tmp_jet_index_candidate.append(jet[k][1])
                            tmp_jet_index_candidate.append(b_jet[j][1])
                            tmp_jet_index_candidate.append(jet[k][2])
                            tmp_jet_index_candidate.append(jet[k][3])
                            tmp_jet_index_candidate.append(b_jet[j][2])
                            tmp_jet_index_candidate.append(b_jet[j][3])
                            _jet_index_candidate.append(np.array(tmp_jet_index_candidate))
                    _jet_index_candidate = np.array(_jet_index_candidate)
            elif lack_of_bjet < 0:
                b_jet = []
                jet = []
                _bjet = itertools.combinations(tmp_bjet_list, 4)
        
                for a in _bjet:
                    b_jet.append(a)

                b_jet = np.array(b_jet, dtype='object')

                _jet_index_candidate = []

                for i in range(len(b_jet)):

                    jet = []
                    
                    tmp_jet_index = tmp_bjet_list.copy()
                    for c in range(len(b_jet[i])): 
                        
                        _tmp = b_jet[i][c]
                        tmp_jet_index.remove(_tmp)
                    
                    _tmp_jet_list = tmp_jet_list.copy()

                    for d in tmp_jet_index:
                        _tmp_jet_list.append(d)
                    
                    _jet = itertools.permutations(_tmp_jet_list, 4)

                    for b in _jet:
                        jet.append(b)

                    jet = np.array(jet, dtype='object')
                    
                    for j in range(len(jet)):
                    
                        tmp_jet_index_candidate = []
                        tmp_jet_index_candidate.append(b_jet[i][0])
                        tmp_jet_index_candidate.append(jet[j][0])
                        tmp_jet_index_candidate.append(jet[j][1])
                        tmp_jet_index_candidate.append(b_jet[i][1])
                        tmp_jet_index_candidate.append(jet[j][2])
                        tmp_jet_index_candidate.append(jet[j][3])
                        tmp_jet_index_candidate.append(b_jet[i][2])
                        tmp_jet_index_candidate.append(b_jet[i][3])
                        _jet_index_candidate.append(np.array(tmp_jet_index_candidate))
                _jet_index_candidate = np.array(_jet_index_candidate)
            else: 
                _jet_index_candidate = np.full(8, -1)
    elif MODEL == 'four_top':
        tmp_jet_list = []
        tmp_bjet_list = []

        _jet_index = []

        for i in range(NUM_OF_JETS):
            _jet_index.append(i)

        for i in range(NUM_OF_JETS):
            if BTAG_JET[i] == 1:
                tmp_bjet_list.append(i)
            else :
                tmp_jet_list.append(i)
        if num_of_btag == 4:
            _bjet = itertools.combinations(tmp_bjet_list, 4)
            _jet = itertools.combinations(tmp_jet_list, 8)
            b_jet = np.array([x for x in _bjet])
            jet = np.array([x for x in _jet])
            
            if len(b_jet) < 1e+4 and len(jet) < 1e+4 :
                _jet_index_candidate = []
                for i in range(len(b_jet)):
                    for j in range(len((jet))):
                        tmp_jet_index_candidate = []
                        tmp_jet_index_candidate.append(b_jet[i][0])
                        tmp_jet_index_candidate.append(jet[j][0])
                        tmp_jet_index_candidate.append(jet[j][1])
                        tmp_jet_index_candidate.append(b_jet[i][1])
                        tmp_jet_index_candidate.append(jet[j][2])
                        tmp_jet_index_candidate.append(jet[j][3])
                        tmp_jet_index_candidate.append(b_jet[i][2])
                        tmp_jet_index_candidate.append(jet[j][4])
                        tmp_jet_index_candidate.append(jet[j][5])
                        tmp_jet_index_candidate.append(b_jet[i][3])
                        tmp_jet_index_candidate.append(jet[j][6])
                        tmp_jet_index_candidate.append(jet[j][7])
                        _jet_index_candidate.append(np.array(tmp_jet_index_candidate))
                _jet_index_candidate = np.array(_jet_index_candidate)
            else: 
                _jet_index_candidate = np.full(12, -1)
        elif num_of_btag != 4:
            
            require_num_bjet = 4
            lack_of_bjet = require_num_bjet - num_of_btag 
            if lack_of_bjet > 0:
                _jet_for_append = itertools.combinations(tmp_jet_list, lack_of_bjet)
                
                jet_for_append = np.array([x for x in _jet_for_append], dtype='object')
                
                for i in range(len(jet_for_append)):

                    tmp_jet_index = tmp_jet_list.copy()
                    for c in range(len(jet_for_append[i])): 
                        _tmp = jet_for_append[i][c]
                        tmp_jet_index.remove(_tmp)

                    tmp_bjet_index = tmp_bjet_list.copy()

                    for d in jet_for_append[i]:                    
                        tmp_bjet_index.append(int(d))

                    _bjet = itertools.permutations(tmp_bjet_index, 4)
                    _jet = itertools.permutations(tmp_jet_index, 8)

                    jet = np.array([x for x in _jet])
                    b_jet = np.array([x for x in _bjet])

                    if len(b_jet) < 1e+4 and len(jet) < 1e+4:
                        _jet_index_candidate = []
                        for j in range(len(b_jet)):
                            for k in range(len((jet))):
                                tmp_jet_index_candidate = []
                                tmp_jet_index_candidate.append(b_jet[j][0])
                                tmp_jet_index_candidate.append(jet[k][0])
                                tmp_jet_index_candidate.append(jet[k][1])
                                tmp_jet_index_candidate.append(b_jet[j][1])
                                tmp_jet_index_candidate.append(jet[k][2])
                                tmp_jet_index_candidate.append(jet[k][3])
                                tmp_jet_index_candidate.append(b_jet[j][2])
                                tmp_jet_index_candidate.append(jet[k][4])
                                tmp_jet_index_candidate.append(jet[k][5])
                                tmp_jet_index_candidate.append(b_jet[j][3])
                                tmp_jet_index_candidate.append(jet[k][6])
                                tmp_jet_index_candidate.append(jet[k][7])
                                _jet_index_candidate.append(np.array(tmp_jet_index_candidate))
                        _jet_index_candidate = np.array(_jet_index_candidate)
                    else: 
                        _jet_index_candidate = np.full(12, -1)
            elif lack_of_bjet < 0:
                _bjet = itertools.combinations(tmp_bjet_list, 4)

                b_jet = np.array([x for x in _bjet], dtype='object')
                
                if len(b_jet) < 1e+4:
                    _jet_index_candidate = []
                    for i in range(len(b_jet)):

                        tmp_jet_index = tmp_bjet_list.copy()
                        for c in range(len(b_jet[i])): 

                            _tmp = b_jet[i][c]
                            tmp_jet_index.remove(_tmp)

                        _tmp_jet_list = tmp_jet_list.copy()

                        for d in tmp_jet_index:
                            _tmp_jet_list.append(d)
                        
                        _jet = itertools.permutations(_tmp_jet_list, 8)

                        jet = np.array([x for x in _jet], dtype='object')
                        if len(jet) < 1e+4:
                            for j in range(len(jet)):
                                tmp_jet_index_candidate = []
                                tmp_jet_index_candidate.append(b_jet[i][0])
                                tmp_jet_index_candidate.append(jet[j][0])
                                tmp_jet_index_candidate.append(jet[j][1])
                                tmp_jet_index_candidate.append(b_jet[i][1])
                                tmp_jet_index_candidate.append(jet[j][2])
                                tmp_jet_index_candidate.append(jet[j][3])
                                tmp_jet_index_candidate.append(b_jet[i][2])
                                tmp_jet_index_candidate.append(jet[j][4])
                                tmp_jet_index_candidate.append(jet[j][5])
                                tmp_jet_index_candidate.append(b_jet[i][3])
                                tmp_jet_index_candidate.append(jet[j][6])
                                tmp_jet_index_candidate.append(jet[j][7])
                                _jet_index_candidate.append(np.array(tmp_jet_index_candidate))
                            else: 
                                tmp_jet_index_candidate = np.full(12, -1)
                                _jet_index_candidate.append(np.array(tmp_jet_index_candidate))
                    _jet_index_candidate = np.array(_jet_index_candidate)
                else: 
                    _jet_index_candidate = np.full(12, -1)
        else:
            pass
    else:
        print("Please select a correct model.")
    return _jet_index_candidate

@nb.jit(nopython=True)
def calc_chi2(CAND_LIST: np.ndarray, PT: np.ndarray, ETA: np.ndarray, PHI: np.ndarray, MASS: np.ndarray, EXTRA: str, MODEL: str) -> Union[float, np.ndarray]:
    pointer = np.array([-1.0, 1.0, 1.0, 1.0])
    if MODEL == "ttbar":
        min_chi2 = -1
        m_W = 81.3
        m_top = 172.7
        sigma_W = 12.3
        sigma_t = 26.3
        set_kinematics = []
        for i in range(0, len(CAND_LIST), 3):
            subset = CAND_LIST[i: i+3]
            _kinematics = []
            for a in subset:
                _px = PT[a]*np.cos(PHI[a])
                _py = PT[a]*np.sin(PHI[a])
                _pz = PT[a]*np.sinh(ETA[a])
                _mass = MASS[a]
                _e = np.sqrt((_px**2 + _py**2 + _pz**2) + _mass**2)
                _kinematics.append([_px, _py, _pz, _e])
            set_kinematics.append(_kinematics)
        W_plus_inv = np.sqrt(np.sum(((np.array(set_kinematics[0][1]) + np.array(set_kinematics[0][2]))*pointer)))
        W_minus_inv = np.sqrt(np.sum(((np.array(set_kinematics[1][1]) + np.array(set_kinematics[1][2]))*pointer)))
        top_inv = np.sqrt(np.sum(((np.array(set_kinematics[0][0]) + np.array(set_kinematics[0][1]) + np.array(set_kinematics[0][2]))*pointer)))
        anti_top_inv = np.sqrt(np.sum(((np.array(set_kinematics[1][0]) + np.array(set_kinematics[1][1]) + np.array(set_kinematics[1][2]))*pointer)))
        if EXTRA == "normal":
            chi2_part_1 = (top_inv - anti_top_inv)**2
            chi2_part_2 = (W_plus_inv - m_W)**2
            chi2_part_3 = (W_minus_inv - m_W)**2
            chi2_tmp = chi2_part_1/(2*(sigma_t**2)) + chi2_part_2/sigma_W**2 + chi2_part_3/sigma_W**2
        elif EXTRA == "mtop_base":
            chi2_part_1 = (top_inv - m_top)**2
            chi2_part_2 = (anti_top_inv - m_top)**2
            chi2_part_3 = (W_plus_inv - m_W)**2
            chi2_part_4 = (W_minus_inv - m_W)**2
            chi2_tmp = chi2_part_1/(sigma_t**2) + chi2_part_2/(sigma_t**2)  + chi2_part_3/(sigma_W**2) + chi2_part_4/(sigma_W**2)
        else: 
            print("Please input a available extra option")
    elif MODEL == 'ttH':
        min_chi2 = -1
        m_W = 81.3
        m_h = 125
        sigma_W = 18.7
        sigma_t = 28.8
        sigma_h = 22.3
        set_kinematics = []
        for i in range(0, len(CAND_LIST), 3):
            subset = CAND_LIST[i: i+3]
            _kinematics = []
            for a in subset:
                _px = PT[a]*np.cos(PHI[a])
                _py = PT[a]*np.sin(PHI[a])
                _pz = PT[a]*np.sinh(ETA[a])
                _mass = MASS[a]
                _e = np.sqrt((_px**2 + _py**2 + _pz**2) + _mass**2)
                _kinematics.append([_px, _py, _pz, _e])
            set_kinematics.append(_kinematics)
        _kinematics = []
        for a in  CAND_LIST[6: 9]:
            _px = PT[a]*np.cos(PHI[a])
            _py = PT[a]*np.sin(PHI[a])
            _pz = PT[a]*np.sinh(ETA[a])
            _mass = MASS[a]
            _e = np.sqrt((_px**2 + _py**2 + _pz**2) + _mass**2)
            _kinematics.append([_px, _py, _pz, _e])
        set_kinematics.append(_kinematics)
        W_plus_inv = np.sqrt(np.sum(((np.array(set_kinematics[0][1]) + np.array(set_kinematics[0][2]))*pointer)))
        W_minus_inv = np.sqrt(np.sum(((np.array(set_kinematics[1][1]) + np.array(set_kinematics[1][2]))*pointer)))
        top_inv = np.sqrt(np.sum(((np.array(set_kinematics[0][0]) + np.array(set_kinematics[0][1]) + np.array(set_kinematics[0][2]))*pointer)))
        anti_top_inv = np.sqrt(np.sum(((np.array(set_kinematics[1][0]) + np.array(set_kinematics[1][1]) + np.array(set_kinematics[1][2]))*pointer)))
        higgs_inv = np.sqrt(np.sum(((np.array(set_kinematics[2][0]) + np.array(set_kinematics[2][1]))*pointer)))
        chi2_part_1 = (top_inv - anti_top_inv)**2
        chi2_part_2 = (W_plus_inv - m_W)**2
        chi2_part_3 = (W_minus_inv - m_W)**2
        chi2_part_4 = (higgs_inv - m_h)**2
        chi2_tmp = chi2_part_1/(2*(sigma_t**2)) + chi2_part_2/sigma_W**2 + chi2_part_3/sigma_W**2 + (chi2_part_4)/sigma_h**2
    elif MODEL == 'four_top':
        min_chi2 = -1
        m_W = 81.3
        m_top = 172.7
        sigma_W = 18.7
        sigma_t = 28.8
        set_kinematics = []
        for i in range(0, len(CAND_LIST), 3):
            subset = CAND_LIST[i: i+3]
            _kinematics = []
            for a in subset:
                _py = PT[a]*np.sin(PHI[a])
                _pz = PT[a]*np.sinh(ETA[a])
                _mass = MASS[a]
                _e = np.sqrt((_px**2 + _py**2 + _pz**2) + _mass**2)
                _kinematics.append([_px, _py, _pz, _e])
            set_kinematics.append(_kinematics)
        W_plus_inv_1 = np.sqrt(np.sum(((np.array(set_kinematics[0][1]) + np.array(set_kinematics[0][2]))*pointer)))
        W_minus_inv_1 = np.sqrt(np.sum(((np.array(set_kinematics[1][1]) + np.array(set_kinematics[1][2]))*pointer)))
        W_plus_inv_2 = np.sqrt(np.sum(((np.array(set_kinematics[2][1]) + np.array(set_kinematics[2][2]))*pointer)))
        W_minus_inv_2 = np.sqrt(np.sum(((np.array(set_kinematics[3][1]) + np.array(set_kinematics[3][2]))*pointer)))
        top_inv_1 = np.sqrt(np.sum(((np.array(set_kinematics[0][0]) + np.array(set_kinematics[0][1]) + np.array(set_kinematics[0][2]))*pointer)))
        anti_top_inv_1 = np.sqrt(np.sum(((np.array(set_kinematics[1][0]) + np.array(set_kinematics[1][1]) + np.array(set_kinematics[1][2]))*pointer)))
        top_inv_2 = np.sqrt(np.sum(((np.array(set_kinematics[2][0]) + np.array(set_kinematics[2][1]) + np.array(set_kinematics[2][2]))*pointer)))
        anti_top_inv_2 = np.sqrt(np.sum(((np.array(set_kinematics[3][0]) + np.array(set_kinematics[3][1]) + np.array(set_kinematics[3][2]))*pointer)))
        chi2_part_1 = (top_inv_1 - m_top)**2
        chi2_part_2 = (top_inv_2 - m_top)**2
        chi2_part_3 = (W_plus_inv_1 - m_W)**2
        chi2_part_4 = (W_minus_inv_1 - m_W)**2
        chi2_part_5 = (anti_top_inv_1 - m_top)**2
        chi2_part_6 = (anti_top_inv_2 - m_top)**2
        chi2_part_7 = (W_plus_inv_2 - m_W)**2
        chi2_part_8 = (W_minus_inv_2 - m_W)**2
        chi2_tmp =  chi2_part_1/sigma_t**2 + chi2_part_2/sigma_t**2 \
                    + chi2_part_3/sigma_W**2 + chi2_part_4/sigma_W**2 \
                    + chi2_part_5/sigma_t**2 + chi2_part_6/sigma_t**2 \
                    + chi2_part_7/sigma_W**2 + chi2_part_8/sigma_W**2
    else: 
        print("Please select a correct model.")

    return chi2_tmp, CAND_LIST

def chi_square_minimizer(jet_pt_chi2: np.ndarray, 
                         jet_eta_chi2: np.ndarray, 
                         jet_phi_chi2: np.ndarray, 
                         jet_btag_chi2: np.ndarray, 
                         jet_mass_chi2: np.ndarray, 
                         MODEL: str, 
                         EXTRA: str, 
                         NUM_OF_JETS: int, 
                        )-> Union[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    This is a function for doing chi-square reconstruction.
    PARTON_ETA: Array, a list of partons's eta in a event.
    PARTON_PHI: Array, a list of partons's phi in a event.
    JET_ETA: Array, a list of jet's eta in a event.
    JET_PHI: Array, a list of jet's phi in a event.
    """
    jet_btag_chi2 = jet_btag_chi2[:NUM_OF_JETS]
    num_of_btag = np.sum(np.array(jet_btag_chi2) ==1)
    
    if MODEL == 'ttbar':
        _parton_jet_index = np.full(6, -1)       
        jet_index_candidate = gen_combinaton(NUM_OF_JETS, MODEL, jet_btag_chi2)
        _cand_record  = []
        _chi2_value = []
        for i in range(len(jet_index_candidate)):
            tmp_chi2_value, tmp_cand = calc_chi2(jet_index_candidate[i], jet_pt_chi2, jet_eta_chi2, jet_phi_chi2, jet_mass_chi2, EXTRA, MODEL)         
            _cand_record.append(tmp_cand)
            _chi2_value.append(tmp_chi2_value)
        _chi2_value = np.array(_chi2_value)
        _cand_record = np.array([x for x in _cand_record])

        not_minus_1 = _chi2_value != -1

        _chi2_value =  _chi2_value[not_minus_1]
        smallest_10_idx = np.argsort(_chi2_value)[:10]

        smallest_10_chi2_value = _chi2_value[smallest_10_idx]
        if len(smallest_10_chi2_value) < 10:
            origin_len = len(smallest_10_chi2_value)
            smallest_10_chi2_candidate = np.array([ x for x in _cand_record[smallest_10_idx]], dtype=np.int8)
            smallest_10_chi2_value = np.pad(smallest_10_chi2_value, (0, 10-origin_len), 'constant', constant_values=(0, -1) )
            smallest_10_chi2_candidate = np.pad(smallest_10_chi2_candidate, ((0, 10-origin_len), (0, 0)), 'constant', constant_values=(0, -1) )
        else:
            smallest_10_chi2_candidate = np.array([ x for x in _cand_record[smallest_10_idx]], dtype=np.int8)
        for a in range(len(smallest_10_chi2_value)):
            if np.isnan(smallest_10_chi2_value[a]) == True or smallest_10_chi2_value[a] == 'nan':
                smallest_10_chi2_value[a] = -1
        min_chi2 = smallest_10_chi2_value[0]

        _jet_parton_index = np.full(len(jet_pt_chi2), -1, dtype=np.int8)
        _parton_jet_index = np.array(_cand_record[smallest_10_idx[0]], dtype=np.int8)

        for k in range(NUM_OF_JETS):
            for l in range(len(_parton_jet_index)):
                if _parton_jet_index[l] == int(k):
                    _jet_parton_index[k] = int(l)
                else :
                    pass
        del _chi2_value, _cand_record 
        return min_chi2, _parton_jet_index, _jet_parton_index, smallest_10_chi2_candidate, smallest_10_chi2_value
    elif MODEL == 'ttH':
        _parton_jet_index = np.full(8, -1)
        jet_index_candidate = gen_combinaton(NUM_OF_JETS, MODEL, jet_btag_chi2)
        _cand_record  = []
        _chi2_value = []
        for i in range(len(jet_index_candidate)):
            tmp_chi2_value, tmp_cand = calc_chi2(jet_index_candidate[i], jet_pt_chi2, jet_eta_chi2, jet_phi_chi2, jet_mass_chi2, EXTRA, MODEL)         
            _cand_record.append(tmp_cand)
            _chi2_value.append(tmp_chi2_value)

        _chi2_value = np.array(_chi2_value)
        _cand_record = np.array([x for x in _cand_record])
#         print(_chi2_value, _cand_record)
        not_minus_1 = _chi2_value != -1
        _chi2_value =  _chi2_value[not_minus_1]
        smallest_10_idx = np.argsort(_chi2_value)[:10]

        smallest_10_chi2_value = _chi2_value[smallest_10_idx]
        if len(smallest_10_chi2_value) < 10:
            origin_len = len(smallest_10_chi2_value)
            smallest_10_chi2_candidate = np.array([ x for x in _cand_record[smallest_10_idx]], dtype=np.int8)
            smallest_10_chi2_value = np.pad(smallest_10_chi2_value, (0, 10-origin_len), 'constant', constant_values=(0, -1) )
            smallest_10_chi2_candidate = np.pad(smallest_10_chi2_candidate, ((0, 10-origin_len), (0, 0)), 'constant', constant_values=(0, -1) )
        else:
            smallest_10_chi2_candidate = np.array([ x for x in _cand_record[smallest_10_idx]], dtype=np.int8)
        for a in range(len(smallest_10_chi2_value)):
            if np.isnan(smallest_10_chi2_value[a]) == True or smallest_10_chi2_value[a] == 'nan':
                smallest_10_chi2_value[a] = -1
        min_chi2 = smallest_10_chi2_value[0]
        
        _jet_parton_index = np.full(len(jet_pt_chi2), -1, dtype=np.int8)
        _parton_jet_index = np.array(_cand_record[smallest_10_idx[0]], dtype=np.int8)

        for k in range(NUM_OF_JETS):
            for l in range(len(_parton_jet_index)):
                if _parton_jet_index[l] == int(k):
                    _jet_parton_index[k] = int(l)
                else :
                    pass

        del _chi2_value, _cand_record 
        return min_chi2, _parton_jet_index, _jet_parton_index, smallest_10_chi2_candidate, smallest_10_chi2_value
    elif MODEL == "four_top":
        _parton_jet_index = np.full(12, -1)
        jet_index_candidate = gen_combinaton(NUM_OF_JETS, MODEL, jet_btag_chi2)
        if len(jet_index_candidate) < 1e+5 and (np.array(jet_index_candidate) == -1).sum() == 0 and len(jet_index_candidate) >= 1:
            _cand_record  = []
            _chi2_value = []
            for i in tqdm.trange(len(jet_index_candidate)):
                tmp_chi2_value, tmp_cand = calc_chi2(jet_index_candidate[i], jet_pt_chi2, jet_eta_chi2, jet_phi_chi2, jet_mass_chi2, EXTRA, MODEL)         
                _cand_record.append(tmp_cand)
                _chi2_value.append(tmp_chi2_value)

            _chi2_value = np.array(_chi2_value)
            _cand_record = np.array([x for x in _cand_record])
            not_minus_1 = _chi2_value != -1

            _chi2_value =  _chi2_value[not_minus_1]
            smallest_10_idx = np.argsort(_chi2_value)[:10]

            smallest_10_chi2_value = _chi2_value[smallest_10_idx]
            if len(smallest_10_chi2_value) < 10:
                origin_len = len(smallest_10_chi2_value)
                smallest_10_chi2_candidate = np.array([ x for x in _cand_record[smallest_10_idx]], dtype=np.int8)
                smallest_10_chi2_value = np.pad(smallest_10_chi2_value, (0, 10-origin_len), 'constant', constant_values=(0, -1) )
                smallest_10_chi2_candidate = np.pad(smallest_10_chi2_candidate, ((0, 10-origin_len), (0, 0)), 'constant', constant_values=(0, -1) )
            else:
                smallest_10_chi2_candidate = np.array([ x for x in _cand_record[smallest_10_idx]], dtype=np.int8)
            for a in range(len(smallest_10_chi2_value)):
                if np.isnan(smallest_10_chi2_value[a]) == True or smallest_10_chi2_value[a] == 'nan':
                    smallest_10_chi2_value[a] = -1
            min_chi2 = smallest_10_chi2_value[0]

            _jet_parton_index = np.full(len(jet_pt_chi2), -1, dtype=np.int8)
            _parton_jet_index = np.array(_cand_record[smallest_10_idx[0]], dtype=np.int8)

            for k in range(NUM_OF_JETS):
                for l in range(len(_parton_jet_index)):
                    if _parton_jet_index[l] == int(k):
                        _jet_parton_index[k] = int(l)
                    else :
                        pass

            del _chi2_value, _cand_record 
        else:
            min_chi2 = -1
            _parton_jet_index = np.full(12, -1)
            _jet_parton_index = np.full(12, -1)
            smallest_10_chi2_value = np.full(10, -1) 
            smallest_10_chi2_candidate = np.full((10, 12), -1)
        return min_chi2, _parton_jet_index, _jet_parton_index, smallest_10_chi2_candidate, smallest_10_chi2_value

# def purity_classifier(prediction, truth_match, mode, model):
#     if model == 'ttbar':
#         if mode == "pair":
#             left_truth_match = truth_match[:3]
#             right_truth_match = truth_match[3:]

#             left_prediction = prediction[:3]    
#             right_prediction = prediction[3:]    

#             b_1_truth_match = left_truth_match[0]
#             b_2_truth_match = right_truth_match[0]
#             j_12_truth_match = set(left_truth_match[1:])
#             j_34_truth_match = set(right_truth_match[1:])
#             truth_match_b_pair = [b_1_truth_match, b_2_truth_match]

#             b_1_prediction = left_prediction[0]
#             b_2_prediction = right_prediction[0]
#             j_12_prediction = set(left_prediction[1:])
#             j_34_prediction = set(right_prediction[1:])
#             prediction_b_pair = [b_1_prediction, b_2_prediction]
            
#             if truth_match_b_pair[0] == prediction_b_pair[0] and truth_match_b_pair[1] == prediction_b_pair[1]:
#                 _case = 1
#             elif truth_match_b_pair[0] == prediction_b_pair[1] and truth_match_b_pair[1] == prediction_b_pair[0]:
#                 _case = 2
#             elif truth_match_b_pair[0] == prediction_b_pair[0] and truth_match_b_pair[1] != prediction_b_pair[1]:
#                 _case = 3
#             elif truth_match_b_pair[0] != prediction_b_pair[0] and truth_match_b_pair[1] == prediction_b_pair[1]:
#                 _case = 4
#             elif truth_match_b_pair[0] == prediction_b_pair[1] and truth_match_b_pair[1] != prediction_b_pair[0]:
#                 _case = 5
#             elif truth_match_b_pair[0] != prediction_b_pair[1] and truth_match_b_pair[1] == prediction_b_pair[0]:
#                 _case = 6
#             else: 
#                 _case = 7
            
#             if _case == 1:
#                 if j_12_truth_match == j_12_prediction and j_34_truth_match == j_34_prediction:
#                     _correct_left = _correct_right = 1
#                 elif j_12_truth_match != j_12_prediction and j_34_truth_match == j_34_prediction:
#                     _correct_left = 0 
#                     _correct_right = 1
#                 elif j_12_truth_match == j_12_prediction and j_34_truth_match != j_34_prediction:
#                     _correct_left = 1
#                     _correct_right = 0
#                 elif j_12_truth_match != j_12_prediction and j_34_truth_match != j_34_prediction:
#                     _correct_left = _correct_right = 0
#                 else:
#                     print("Error occur in function.")

#             elif _case == 2:
#                 if j_12_truth_match == j_34_prediction and j_34_truth_match == j_12_prediction:
#                     _correct_left = _correct_right = 1
#                 elif j_12_truth_match != j_34_prediction and j_34_truth_match == j_12_prediction:
#                     _correct_left = 0 
#                     _correct_right = 1
#                 elif j_12_truth_match == j_34_prediction and j_34_truth_match != j_12_prediction:
#                     _correct_left = 1
#                     _correct_right = 0
#                 elif j_12_truth_match != j_34_prediction and j_34_truth_match != j_12_prediction:
#                     _correct_left = _correct_right = 0
#                 else:
#                     print("Error occur in function.")
#             elif _case == 3:
#                 _correct_right = 0
#                 if j_12_truth_match == j_12_prediction:
#                     _correct_left = 1
#                 else:
#                     _correct_left = 0
#             elif _case == 4:
#                 _correct_left = 0
#                 if j_34_truth_match == j_34_prediction:
#                     _correct_right = 1
#                 else:
#                     _correct_right = 0
#             elif _case == 5:
#                 _correct_right = 0
#                 if j_12_truth_match == j_34_prediction:
#                     _correct_left = 1
#                 else:
#                     _correct_left = 0
#             elif _case == 6:
#                 _correct_left = 0
#                 if j_34_truth_match == j_12_prediction:
#                     _correct_right = 1
#                 else:
#                     _correct_right = 0
#             elif _case == 7:
#                 _correct_left  = _correct_right = 0
#             else: 
#                 pass

#             return _correct_left, _correct_right
        
#         elif mode == "left":
#             left_truth_match = set(truth_match[:3])
#             left_prediction = set(prediction[:3])
#             right_prediction = set(prediction[3:])

#             if left_truth_match == left_prediction: 
#                 _correct_left = 1
#             elif left_truth_match == right_prediction: 
#                 _correct_left = 1
#             else: 
#                 _correct_left = 0
                
#             return _correct_left
#         elif mode == "right":
#             right_truth_match = set(truth_match[3:])
#             right_prediction = set(prediction[3:])
#             left_prediction = set(prediction[:3])
            
#             if right_truth_match == right_prediction: 
#                 _correct_right = 1
#             elif  right_truth_match == left_prediction: 
#                 _correct_right = 1
#             else: 
#                 _correct_right = 0
#             return  _correct_right
#         else: 
#             print("Error occur in function")
#     else: 
#         print("Please select a available model. (Only ttbar model(fully hadronic decay) available currently.)")

