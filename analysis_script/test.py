#Import package
import os, sys, time
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
                            "Particle.Mass"]) 
        self.num_of_event = len(self.event)
        self.jet = root2array(self.file, "Delphes;1", branches=["Jet.PT", "Jet.Eta", "Jet.Phi","Jet.Mass"])
        self._Status, self._M1, self._M2, self._D1, self._D2, self._PID, self._PT, self._Eta, self._Phi, self._Mass = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 
        self.Labels = ["Status", "M1", "M2", "D1", "D2", "PID", "PT", "Eta", "Phi", "Mass"]
        
    def parse_event(self):

        return self.event

    def parse_jet(self):
        
        return self.jet

    def count_particles(self):           

        self.new_array = np.unique( self.event[0][self._PID] ) 
        self.count = len(self.new_array)
        return self.count

    def in_and_out_finder(self):

        print("Work in progress")
    

    def save_parsed_file(self):
        for i in range(len(self.num_of_event)):
            time_tag =  time.strftime("%Y_%m_%d_%H_%M", time.localtime())
            FILENAME = ""
            seq = ("parsed_", time_tag, "_event_", i, ".txt")
            os.system("mkdir {0}/parsed_event_{1}".format(self.path, time_tag))
            store_path = os.join.path(self.path, "parsed_event_{0}".format(time_tag),FILENAME.join(seq))

            idx = np.linspace(0, len( self.event[0][0])-1, num = len( self.event[0][0]) )
            patron_dict = {
                "Index": idx,
                "Status":  self.event[i][0],
                "Mother_1":  self.event[i][1],
                "Mother_2":  self.event[i][2],
                "Daughter_1":  self.event[i][3],
                "Daughter_2":  self.event[i][4],
                "PID":  self.event[i][5],
                "PT":  self.event[i][6],
                "Eta":  self.event[i][7],
                "Phi":  self.event[i][8],
                "Mass":  self.event[i][9]
                }
            patron_df = pd.DataFrame(patron_dict)
            patron_df.to_csv(path_or_buf=store_path, header=True)
            

            # with open(store_path, "w") as f:
            #     f.write("Index,status,M1,M2,D1,D2,PID,PT,ETA,PHI,MASS")
            #     f.write("\n")
            #     for j in range(len(self.event[0][0])):
            #         f.write("{0},{1},{2},{3},{4},{5},{6},{7:.3f},{8:.3f},{9:.3f},{10:.3f}".format( 
            #                 j, 
            #                 self.event[i][self._Status][j], 
            #                 self.event[i][self._M1][j],
            #                 self.event[i][self._M2][j], 
            #                 self.event[i][self._D1][j],  
            #                 self.event[i][self._D2][j], 
            #                 self.event[i][self._PID][j], 
            #                 self.event[i][self._PT][j], 
            #                 self.event[i][self._Eta][j], 
            #                 self.event[i][self._Phi][j],
            #                 self.event[i][self._Mass][j]))
            #         f.write("\n")
            return store_path

def main():

    root_file = sys.argv[1]  #your root file path
    time_tag = time.strftime("%Y_%m_%d_%H_%M", time.localtime())
    
    start = time.time()
    PATH = os.getcwd()
    parsed_event = parse_function(root_file, PATH)

    event_log = parsed_event.parse_event()
    end = time.time() - start
    jet_log = parsed_event.parse_jet()
    particle_numbers = parsed_event.count_particles()
    
    
    print("Cost {0}s.".format(end))
    
    # print(len(idx), len(event_log[0][0]), len(event_log[0][1]), 
    #     len(event_log[0][2]), len(event_log[0][3]), len(event_log[0][4]),
    #     len(event_log[0][5]), len(event_log[0][6]), len(event_log[0][7]),
    #     len(event_log[0][8]), len(event_log[0][9])
    # )

    # patron_dict = {
    #     "Index": idx,
    #     "Status": event_log[0][0],
    #     "Mother_1": event_log[0][1],
    #     "Mother_2": event_log[0][2],
    #     "Daughter_1": event_log[0][3],
    #     "Daughter_2": event_log[0][4],
    #     "PID": event_log[0][5],
    #     "PT": event_log[0][6],
    #     "Eta": event_log[0][7],
    #     "Phi": event_log[0][8],
    #     "Mass": event_log[0][9]
    # }
    
    # patron_df = pd.DataFrame(patron_dict)
    # print(patron_df)

    #print(jet_log)

    #saved_path = parsed_event.save_parsed_file()
    #print("The data has been parsed! There are {0} particles in this events.\nThe result has been store in '{1}'.".format(particle_numbers, saved_path))


if __name__ == "__main__": 
    main()

