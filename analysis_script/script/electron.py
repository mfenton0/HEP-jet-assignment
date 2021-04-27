"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
import uproot
import pandas as pd
import numpy as np
import tqdm

class electron_properties():
    def __init__(self, data, single=True):
        """
        data: str or list, the directory of root file. This variable will become a list. 
        single: boolean, whether loading single file or not.
        """
        if single:
            print("Loading electron information.")
            self.pt = data.array("Electron.PT")
            self.eta = data.array('Electron.Eta')
            self.phi = data.array('Electron.Phi')
        else:
            print("Loading electron information.")
            num_of_files = len(data)
            count = 1
            _data = []
            for a in data:
                try:
                    print(f"Padding root file from {a}. Progress: {count}/{num_of_files}.")
                    tmp = uproot.open(a)['Delphes']
                    _data.append(tmp)
                    count += 1
                except:
                    print('Please check input file path.')
            print("Loading electron information.")
            for i in tqdm.trange(num_of_files):
                if i == 0 :
                    _electron_pt = _data[i].array('Electron.PT')
                    _electron_eta = _data[i].array('Electron.Eta')
                    _electron_phi = _data[i].array('Electron.Phi')
                else: 
                    _electron_pt = np.concatenate((_electron_pt, _data[i].array('Electron.PT')))
                    _electron_eta = np.concatenate((_electron_eta, _data[i].array('Electron.Eta')))
                    _electron_phi = np.concatenate((_electron_phi, _data[i].array('Electron.Phi')))
            self.pt = _electron_pt
            self.eta = _electron_eta
            self.phi = _electron_phi
    def dataframelize(self, index):

        idx = np.linspace(0, len( self.pt[index])-1, num = len( self.pt[index]) )

        electron_dict = {
                "Index": idx,
                "Electron_pt":  self.pt[index],
                "Electron_eta":  self.eta[index],
                "Electron_phi":  self.phi[index],

            }
        electron_df = pd.DataFrame(electron_dict)
        return electron_df
