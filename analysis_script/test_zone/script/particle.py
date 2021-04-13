"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
import uproot
import pandas as pd
import numpy as np
import tqdm

class particle_properties():
    def __init__(self, data, single=True):
        
        """
        data: str or list, the directory of root file. This variable will become a list. 
        single: boolean, whether loading single file or not.
        """
        if single:
            print("Loading particle information.")
            self.event = data.array('Event')
            self.pt = data.array('Particle.PT')
            self.eta = data.array('Particle.Eta')
            self.phi = data.array('Particle.Phi')
            self.pid = data.array('Particle.PID')
            self.M1 = data.array('Particle.M1')
            self.M2 = data.array('Particle.M2')
            self.D1 = data.array('Particle.D1')
            self.D2 = data.array('Particle.D2')
            self.status = data.array('Particle.Status')
            self.rapidity = data.array('Particle.Rapidity')
            self.mass = data.array('Particle.Mass')
            self.charge = data.array('Particle.Charge')
        else: 
            num_of_files = len(data)
            _data = []
            count = 0
            for a in data:
                try:
                    print(f"Padding root file from {a}. Progress: {count}/{num_of_files}.")
                    tmp = uproot.open(a)['Delphes']
                    _data.append(tmp)
                    count += 1
                except:
                    print('Please check input file path.')
            print("Loading particle information.")
            for i in tqdm.trange(num_of_files):
                if i == 0 :
                    _particle_event = _data[i].array('Event')
                    _particle_pt = _data[i].array('Particle.PT')
                    _particle_eta = _data[i].array('Particle.Eta')
                    _particle_phi = _data[i].array('Particle.Phi')
                    _particle_pid = _data[i].array('Particle.PID')
                    _particle_M1 = _data[i].array('Particle.M1')
                    _particle_M2 = _data[i].array('Particle.M2')
                    _particle_D1 = _data[i].array('Particle.D1')
                    _particle_D2 = _data[i].array('Particle.D2')
                    _particle_status = _data[i].array('Particle.Status')
                    _particle_rapidity = _data[i].array('Particle.Rapidity')
                    _particle_mass = _data[i].array('Particle.Mass')
                    _particle_charge = _data[i].array('Particle.Charge')
                else: 
                    _particle_event = np.concatenate((_particle_event, _data[i].array('Event')))
                    _particle_pt = np.concatenate((_particle_pt, _data[i].array('Particle.PT')))
                    _particle_eta = np.concatenate((_particle_eta, _data[i].array('Particle.Eta')))
                    _particle_phi = np.concatenate((_particle_phi, _data[i].array('Particle.Phi')))
                    _particle_pid = np.concatenate((_particle_pid,_data[i].array('Particle.PID')))
                    _particle_M1 = np.concatenate((_particle_M1, _data[i].array('Particle.M1')))
                    _particle_M2 = np.concatenate((_particle_M2, _data[i].array('Particle.M2')))
                    _particle_D1 = np.concatenate((_particle_D1, _data[i].array('Particle.D1')))
                    _particle_D2 = np.concatenate((_particle_D2, _data[i].array('Particle.D2')))
                    _particle_status = np.concatenate((_particle_status, _data[i].array('Particle.Status')))
                    _particle_rapidity = np.concatenate((_particle_rapidity, _data[i].array('Particle.Rapidity')))
                    _particle_mass = np.concatenate((_particle_mass, _data[i].array('Particle.Mass')))
                    _particle_charge = np.concatenate((_particle_charge, _data[i].array('Particle.Charge')))

            
            self.event = _particle_event
            self.pt = _particle_pt
            self.eta = _particle_eta
            self.phi = _particle_phi
            self.pid = _particle_pid
            self.M1 = _particle_M1
            self.M2 = _particle_M2
            self.D1 = _particle_D1
            self.D2 = _particle_D2
            self.status = _particle_status
            self.rapidity = _particle_rapidity
            self.mass = _particle_mass
            self.charge = _particle_charge
    def dataframelize(self, index):

        idx = np.linspace(0, len( self.pt[index])-1, num = len( self.pt[index]) )

        patron_dict = {
                "Index": idx,
                "Status":  self.status[index],
                "Mother_1":  self.M1[index],
                "Mother_2":  self.M2[index],
                "Daughter_1":  self.D1[index],
                "Daughter_2":  self.D2[index],
                "PID":  self.pid[index],
                "PT":  self.pt[index],
                "Eta":  self.eta[index],
                "Phi":  self.phi[index],
                "Mass":  self.mass[index]
            }
        patron_df = pd.DataFrame(patron_dict)
        return patron_df

