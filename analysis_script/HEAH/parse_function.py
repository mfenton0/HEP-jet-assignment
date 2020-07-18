#Import package
import os, time, sys
import pandas as pd
from root_numpy import root2array
import numpy as np


class parse_function():

    def __init__(self, to_be_parsed_file, path):

        self.file = to_be_parsed_file
        self.path = path
        #Load root file
        self.event = root2array(self.file, "Delphes;1", 
                branches=["Particle.Status", "Particle.M1",
                            "Particle.M2", "Particle.D1", "Particle.D2","Particle.PID",
                            "Particle.PT","Particle.Eta","Particle.Phi",
                            "Particle.Mass"], start=0, stop=2000, step=None) # in this example I set
        self._Status, self._M1, self._M2, self._D1, self._D2, self._PID, self._PT, self._Eta, self._Phi, self._Mass = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 
        self.Labels = ["Status", "M1", "M2", "D1", "D2", "PID", "PT", "Eta", "Phi", "Mass"]

    def parse_event(self):

        return self.event

    def count_particles(self):           

        self.new_array = np.unique( self.event[0][self._PID] ) 
        self.count = len(self.new_array)
        return self.count

    def check_state(self):

        print("Work in progress")

    def save_parsed_file(self):

        with open(self.path, "w") as f:
            f.write("Index,status,M1,M2,D1,D2,PID,PT,ETA,PHI,MASS")
            f.write("\n")
            for j in range(len(self.event[0][0])):
                f.write("{0},{1},{2},{3},{4},{5},{6:.3f},{7:.3f},{8:.3f},{9:.3f}".format( 
                        j, 
                        self.event[0][self._Status][j], 
                        self.event[0][self._M1][j], 
                        self.event[0][self._D1][j],  
                        self.event[0][self._D2][j], 
                        self.event[0][self._PID][j], 
                        self.event[0][self._PT][j], 
                        self.event[0][self._Eta][j], 
                        self.event[0][self._Phi][j],
                        self.event[0][self._Mass][j]))
                f.write("\n")
        return self.path



