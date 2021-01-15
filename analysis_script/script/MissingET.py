import uproot
import pandas as pd
import numpy as np

class Missing_ET():
    def __init__(self, data):
        self.MET = data.array('MissingET.MET')
        self.eta = data.array('MissingET.Eta')
        self.phi = data.array('MissingET.Phi')
        
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