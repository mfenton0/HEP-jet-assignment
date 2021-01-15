import uproot
import pandas as pd
import numpy as np

class Lepton():
    def __init__(self, data):
        self.electron_pt = data.array("Electron.PT")
        self.electron_eta = data.array('Electron.Eta')
        self.electron_phi = data.array('Electron.Phi')
        self.Muon_pt = data.array("Electron.PT")
        self.Muon_eta = data.array('Muon.Eta')
        self.Muon_phi = data.array('Muon.Phi')
        
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