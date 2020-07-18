#Import package
import os
import time
import sys
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
            # f.write("Index", "\t","Status", "\t","M1"
            #     "\t","M2" ,"\t","D1", "\t","D2", "\t","PID", 
            #     "\t\t","PT" "\t","Eta", "\t\t","Phi", 
            #     "\t\t","Mass")
            f.write("Index,status,M1,M2,D1,D2,PID,PT,ETA,PHI,MASS")
            f.write("\n")
            for j in range(len(self.event[0][0])):
            #for j in range(1000):
                # f.write(j, "\t", event[0][_Status][j],"\t\t",
                #         event[0][_M1][j], "\t", event[0][_M2][j],
                #         "\t", event[0][_D1][j], "\t", event[0][_D2][j],
                #         "\t", str(event[0][_PID][j]).ljust(12, ' '), "\t", round(event[0][_PT][j],0),  "\t",
                #         str(round(event[0][_Eta][j],2)).ljust(12, ' '), "\t",
                #         str(round(event[0][_Phi][j],3)).ljust(12, ' '), "\t",
                #         round(event[0][_Mass][j],3))
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

def main():
    root_file = sys.argv[1]  #your root file path
    time_tag = time.strftime("%Y_%m_%d_%H_%M", time.localtime())
    FILENAME = ""
    seq = ("parsed_", time_tag, ".txt")
    PATH = os.path.join(os.getcwd(), FILENAME.join(seq))
    parsed_event = parse_function(root_file, PATH)
    particle_numbers = parsed_event.count_particles()
    saved_path = parsed_event.save_parsed_file()
    print("The data has been parsed! There are {0} particles in this events. \n The result has been store in '{1}'.".format(particle_numbers, saved_path))


if __name__ == "__main__": 
    main()

