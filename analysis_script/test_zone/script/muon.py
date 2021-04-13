"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
import uproot
import pandas as pd
import numpy as np
import tqdm 

class muon_properties():
    def __init__(self, data, single=True):
        """
        data: str or list, the directory of root file. This variable will become a list. 
        single: boolean, whether loading single file or not.
        """
        if single:
            print("Loading muon information.")
            self.pt = data.array("Muon.PT")
            self.eta = data.array('Muon.Eta')
            self.phi = data.array('Muon.Phi')
        else:
            print("Loading muon information.")
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
            print("Loading muon information.")
            for i in tqdm.trange(num_of_files):
                if i == 0 :
                    _muon_pt = _data[i].array('Muon.PT')
                    _muon_eta = _data[i].array('Muon.Eta')
                    _muon_phi = _data[i].array('Muon.Phi')
                else: 
                    _muon_pt = np.concatenate((_muon_pt, _data[i].array('Muon.PT')))
                    _muon_eta = np.concatenate((_muon_eta, _data[i].array('Muon.Eta')))
                    _muon_phi = np.concatenate((_muon_phi, _data[i].array('Muon.Phi')))
            self.pt = _muon_pt
            self.eta = _muon_eta
            self.phi = _muon_phi
    def dataframelize(self, index):

        idx = np.linspace(0, len( self.pt[index])-1, num = len( self.pt[index]) )

        muon_dict = {
                "Index": idx,
                "Muon_pt":  self.pt[index],
                "Muon_eta":  self.eta[index],
                "Muon_phi":  self.phi[index],
            }
        muon_df = pd.DataFrame(muon_dict)
        return muon_df
