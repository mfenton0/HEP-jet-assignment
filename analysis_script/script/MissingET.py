"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
import uproot
import pandas as pd
import numpy as np
import tqdm

class Missing_ET_properties():
    def __init__(self, data, single=True):
        """
        data: str or list, the directory of root file. This variable will become a list. 
        single: boolean, whether loading single file or not.
        """
        if single:
            print("Loading MissingET information.")
            self.MET = data.array('MissingET.MET')
            self.eta = data.array('MissingET.Eta')
            self.phi = data.array('MissingET.Phi')
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
            print("Loading MissingET information.")
            for i in tqdm.trange(num_of_files):
                if i == 0:
                    _missing_et_met = _data[i].array('MissingET.MET')
                    _missing_et_eta = _data[i].array('MissingET.Eta')
                    _missing_et_phi = _data[i].array('MissingET.Phi')
                else: 
                    _missing_et_met = np.concatenate((_missing_et_met, _data[i].array('MissingET.MET')))
                    _missing_et_eta = np.concatenate((_missing_et_eta, _data[i].array('MissingET.Eta')))
                    _missing_et_phi = np.concatenate((_missing_et_phi, _data[i].array('MissingET.Phi')))
            self.MET = _missing_et_met
            self.eta = _missing_et_eta
            self.phi = _missing_et_phi

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
