import uproot
import numpy as np
import pandas as pd

class jet_properties():
    def __init__(self, data):
        self.event = data.array('Event')
        self.pt = data.array('Jet.PT')
        self.eta = data.array('Jet.Eta')
        self.phi = data.array('Jet.Phi')
        self.btag = data.array('Jet.BTag')
        self.area = data.array('Jet.Area')
        self.mass = data.array('Jet.Mass')
        self.charge = data.array('Jet.Charge')
        self.px = self.pt * np.cos(self.phi)
        self.py = self.pt * np.sin(self.phi)
        self.pz = self.pt * np.sinh(self.eta)
        self.e = np.sqrt( (self.px**2 + self.py**2 + self.pz**2) + self.mass**2 )

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

