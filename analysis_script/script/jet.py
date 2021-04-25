"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
import uproot
import numpy as np
import pandas as pd
import tqdm

class jet_properties():
    def __init__(self, data, single=True):
        """
        data: str or list, the directory of root file. This variable will become a list. 
        single: boolean, whether loading single file or not.
        """
        if single:
            print("Loading jet information.")
            self.event = data.array('Event')
            self.pt = data.array('Jet.PT')
            self.eta = data.array('Jet.Eta')
            self.phi = data.array('Jet.Phi')
            self.btag = data.array('Jet.BTag')
            self.area = data.array('Jet.Area')
            self.mass = data.array('Jet.Mass')
            self.charge = data.array('Jet.Charge')
            self.num_of_jets = data.array('Jet')
        else: 
            num_of_files = len(data)
            count = 0
            _data = []
            for a in data:
                try:
                    print(f"Padding root file from {a}. Progress: {count}/{num_of_files}.")
                    tmp = uproot.open(a)['Delphes']
                    _data.append(tmp)
                    count += 1
                except:
                    print('Please check input file path.')
            print("Loading jet information.")
            for i in tqdm.trange(num_of_files):
                if i == 0 :
                    _jet_event = _data[i].array('Event')
                    _jet_pt = _data[i].array('Jet.PT')
                    _jet_eta = _data[i].array('Jet.Eta')
                    _jet_phi = _data[i].array('Jet.Phi')
                    _jet_btag = _data[i].array('Jet.BTag')
                    _jet_area = _data[i].array('Jet.Area')
                    _jet_mass = _data[i].array('Jet.Mass')
                    _jet_charge = _data[i].array('Jet.Charge')
                    _num_of_jets = _data[i].array('Jet')
                else: 
                    _jet_event = np.concatenate((_jet_event,_data[i].array('Event')))
                    _jet_pt = np.concatenate((_jet_pt, _data[i].array('Jet.PT')))
                    _jet_eta = np.concatenate((_jet_eta, _data[i].array('Jet.Eta')))
                    _jet_phi = np.concatenate((_jet_phi, _data[i].array('Jet.Phi')))
                    _jet_btag = np.concatenate((_jet_btag, _data[i].array('Jet.BTag')))
                    _jet_area = np.concatenate((_jet_area, _data[i].array('Jet.Area')))
                    _jet_mass = np.concatenate((_jet_mass, _data[i].array('Jet.Mass')))
                    _jet_charge = np.concatenate((_jet_charge, _data[i].array('Jet.Charge')))
                    _num_of_jets =  np.concatenate((_num_of_jets, _data[i].array('Jet')))

            self.event = _jet_event
            self.pt = _jet_pt
            self.eta = _jet_eta
            self.phi = _jet_phi
            self.btag = _jet_btag
            self.area = _jet_area
            self.mass = _jet_mass
            self.charge = _jet_charge
            self.num_of_jets = _num_of_jets
    def dataframelize(self, index):

        idx = np.linspace(0, len( self.pt[index])-1, num = len( self.pt[index]) )

        jet_dict = {
                "Index": idx,
                "PT":  self.pt[index],
                "Eta":  self.eta[index],
                "Phi":  self.phi[index],
                "Mass":  self.mass[index],
                "Btag": self.btag[index],
                "Area": self.area[index],
                "Number_of_jet": self.num_of_jets[index],
            }
        jet_df = pd.DataFrame(jet_dict)
        return jet_df

