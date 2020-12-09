"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
import uproot
import pandas as pd 
import numpy as np 
from particle_properties_uproot import particle_properties  #import particle properties helper function from particle_properties.py
from jet_properties_uproot import jet_properties  #import jet properties helper function from jet_properties.py
import h5py, sys, traceback, os, tqdm, configparser
from utilize import delta_R, deltaPhi, pdgid, event_selection, quark_finder, particle_tracing, deltaR_matching, barcode_recorder
import matplotlib.pyplot as plt 

left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height

class properties_for_particle():
    def __init__(self, event, pt, eta, phi, pid, M1, M2, D1, D2, status, rapidity, mass, charge):
        self.event = event
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.pid = pid
        self.M1 = M1
        self.M2 = M2
        self.D1 = D1
        self.D2 = D2
        self.status = status
        self.rapidity = rapidity
        self.mass = mass
        self.charge = charge
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

class properties_for_jet():
    def __init__(self, event, pt, eta, phi, btag, area, mass, charge):
        self.event = event
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.btag = btag
        self.area = area
        self.mass = mass
        self.charge = charge
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

def cutflow(INPUT_FILE, OUTPUT_FILE, MODEL, CONFIG, SINGLE):
    
    config = configparser.ConfigParser()
    config.read(CONFIG)
    cuts = config["Cuts"]
    parameter = config["Parameter"]
    num_of_cuts = int(parameter['NUM_OF_CUTS'])
    pt_cuts = int(parameter['PT_CUTS'])
    eta_cuts = float(parameter['ETA_CUTS'])
    
    if num_of_cuts == 1:
        jet_C1 = int(cuts['cut_1_jet'])
        bjet_C1 = int(cuts['cut_1_bjet'])
    elif num_of_cuts == 2:
        jet_C1 = int(cuts['cut_1_jet'])
        bjet_C1 = int(cuts['cut_1_bjet'])
        jet_C2 = int(cuts['cut_2_jet'])
        bjet_C2 = int(cuts['cut_2_bjet'])
    elif num_of_cuts == 3:
        jet_C1 = int(cuts['cut_1_jet'])
        bjet_C1 = int(cuts['cut_1_bjet'])
        jet_C2 = int(cuts['cut_2_jet'])
        bjet_C2 = int(cuts['cut_2_bjet'])
        jet_C3 = int(cuts['cut_3_jet'])
        bjet_C3 = int(cuts['cut_3_bjet'])
    elif num_of_cuts == 4:
        jet_C1 = int(cuts['cut_1_jet'])
        bjet_C1 = int(cuts['cut_1_bjet'])
        jet_C2 = int(cuts['cut_2_jet'])
        bjet_C2 = int(cuts['cut_2_bjet'])
        jet_C3 = int(cuts['cut_3_jet'])
        bjet_C3 = int(cuts['cut_3_bjet'])
        jet_C4 = int(cuts['cut_4_jet'])
        bjet_C4 = int(cuts['cut_4_bjet'])
    elif num_of_cuts == 5:
        jet_C1 = int(cuts['cut_1_jet'])
        bjet_C1 = int(cuts['cut_1_bjet'])
        jet_C2 = int(cuts['cut_2_jet'])
        bjet_C2 = int(cuts['cut_2_bjet'])
        jet_C3 = int(cuts['cut_3_jet'])
        bjet_C3 = int(cuts['cut_3_bjet'])
        jet_C4 = int(cuts['cut_4_jet'])
        bjet_C4 = int(cuts['cut_4_bjet'])
        jet_C5 = int(cuts['cut_5_jet'])
        bjet_C5 = int(cuts['cut_5_bjet'])
    else:
        print("Please input a number of cuts within 1 ~ 5.")
    

    if int(SINGLE) == 1:
        try:
            print("Loading root file from {0}.".format(INPUT_FILE))
            data  = uproot.open(INPUT_FILE)['Delphes']
        except:
            print('Please check input file path.')

        particle = particle_properties(data)
        jet = jet_properties(data)

        if num_of_cuts == 1:
            marker_event_C1 = []
            marker_jet = []
            marker_bjet = []
            for i in range(len(particle.event)):
                marker_event_C1.append(0)
                marker_jet.append(np.zeros([len(jet.pt[i])]))
                marker_bjet.append(np.zeros([len(jet.pt[i])]))
            
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)

        elif num_of_cuts == 2:
            marker_event_C1 = []
            marker_event_C2 = []
            marker_jet = []
            marker_bjet = []
            for i in range(len(particle.event)):
                marker_event_C1.append(0)
                marker_event_C2.append(0)
                marker_jet.append(np.zeros([len(jet.pt[i])]))
                marker_bjet.append(np.zeros([len(jet.pt[i])]))
            
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_event_C2 = np.asanyarray(marker_event_C2, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)

        elif num_of_cuts == 3:
            marker_event_C1 = []
            marker_event_C2 = []
            marker_event_C3 = []
            marker_jet = []
            marker_bjet = []
            for i in range(len(particle.event)):
                marker_event_C1.append(0)
                marker_event_C2.append(0)
                marker_event_C3.append(0)
                marker_jet.append(np.zeros([len(jet.pt[i])]))
                marker_bjet.append(np.zeros([len(jet.pt[i])]))
            
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_event_C2 = np.asanyarray(marker_event_C2, dtype=object)
            marker_event_C3 = np.asanyarray(marker_event_C3, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)

        elif num_of_cuts == 4:
            marker_event_C1 = []
            marker_event_C2 = []
            marker_event_C3 = []
            marker_event_C4 = []
            marker_jet = []
            marker_bjet = []

            for i in range(len(particle.event)):
                marker_event_C1.append(0)
                marker_event_C2.append(0)
                marker_event_C3.append(0)
                marker_event_C4.append(0)
                marker_jet.append(np.zeros([len(jet.pt[i])]))
                marker_bjet.append(np.zeros([len(jet.pt[i])]))
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_event_C2 = np.asanyarray(marker_event_C2, dtype=object)
            marker_event_C3 = np.asanyarray(marker_event_C3, dtype=object)
            marker_event_C4 = np.asanyarray(marker_event_C4, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)
        elif num_of_cuts == 5:
            marker_event_C1 = []
            marker_event_C2 = []
            marker_event_C3 = []
            marker_event_C4 = []
            marker_event_C5 = []
            marker_jet = []
            marker_bjet = []
            for i in range(len(particle.event)):
                marker_event_C1.append(0)
                marker_event_C2.append(0)
                marker_event_C3.append(0)
                marker_event_C4.append(0)
                marker_event_C5.append(0)
                marker_jet.append(np.zeros([len(jet.pt[i])]))
                marker_bjet.append(np.zeros([len(jet.pt[i])]))
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_event_C2 = np.asanyarray(marker_event_C2, dtype=object)
            marker_event_C3 = np.asanyarray(marker_event_C3, dtype=object)
            marker_event_C4 = np.asanyarray(marker_event_C4, dtype=object)
            marker_event_C5 = np.asanyarray(marker_event_C5, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)
        else:
            print("Please input a number of cuts within 1 ~ 5.")

        print("+------------------------------------------------------------------------------------------------------+")
        print("Applying cuts.")
        print("+------------------------------------------------------------------------------------------------------+")
        for i in tqdm.trange(len(particle.event)):
            for j in range(len(jet.pt[i])):
                if jet.btag[i][j] == 1 and jet.pt[i][j] > pt_cuts and np.abs(jet.eta[i][j]) < eta_cuts:
                    marker_bjet[i][j] = 1 
                else: pass 
            
                if jet.pt[i][j] > pt_cuts and np.abs(jet.eta[i][j]) <= eta_cuts:
                    marker_jet[i][j] = 1
                else: pass 
        print("+------------------------------------------------------------------------------------------------------+")
        print("Cuts applied.")
        print("+------------------------------------------------------------------------------------------------------+")

        print("+------------------------------------------------------------------------------------------------------+")
        print("Computing cutflows.")
        print("+------------------------------------------------------------------------------------------------------+")
        for i in tqdm.trange(len(particle.event)):
            if num_of_cuts == 1:
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
            elif num_of_cuts == 2:
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C2 and np.sum(marker_bjet[i] == 1) >= bjet_C2:
                    marker_event_C2[i] = 1 
            elif num_of_cuts == 3:
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C2 and np.sum(marker_bjet[i] == 1) >= bjet_C2:
                    marker_event_C2[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C3 and np.sum(marker_bjet[i] == 1) >= bjet_C3:
                    marker_event_C3[i] = 1 
            elif num_of_cuts == 4: 
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C2 and np.sum(marker_bjet[i] == 1) >= bjet_C2:
                    marker_event_C2[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C3 and np.sum(marker_bjet[i] == 1) >= bjet_C3:
                    marker_event_C3[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C4 and np.sum(marker_bjet[i] == 1) >= bjet_C4:
                    marker_event_C4[i] = 1 
            elif num_of_cuts == 5:
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C2 and np.sum(marker_bjet[i] == 1) >= bjet_C2:
                    marker_event_C2[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C3 and np.sum(marker_bjet[i] == 1) >= bjet_C3:
                    marker_event_C3[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C4 and np.sum(marker_bjet[i] == 1) >= bjet_C4:
                    marker_event_C4[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C5 and np.sum(marker_bjet[i] == 1) >= bjet_C5:
                    marker_event_C5[i] = 1 
        print("+------------------------------------------------------------------------------------------------------+")
        print("Cutflows computed.")
        print("+------------------------------------------------------------------------------------------------------+")

        print("+------------------------------------------------------------------------------------------------------+")
        print("Start to generate figures.")
        print("+------------------------------------------------------------------------------------------------------+")
        if num_of_cuts == 1:
            cutflow = np.zeros(3)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
        elif num_of_cuts == 2:
            cutflow = np.zeros(4)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
            cutflow[3] = np.sum(marker_event_C2 == 1)
        elif num_of_cuts == 3:
            cutflow = np.zeros(5)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
            cutflow[3] = np.sum(marker_event_C2 == 1)
            cutflow[4] = np.sum(marker_event_C3 == 1)
        elif num_of_cuts == 4:
            cutflow = np.zeros(6)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
            cutflow[3] = np.sum(marker_event_C2 == 1)
            cutflow[4] = np.sum(marker_event_C3 == 1)
            cutflow[5] = np.sum(marker_event_C4 == 1)
        elif num_of_cuts == 5:
            cutflow = np.zeros(7)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
            cutflow[3] = np.sum(marker_event_C2 == 1)
            cutflow[4] = np.sum(marker_event_C3 == 1)
            cutflow[5] = np.sum(marker_event_C4 == 1)
            cutflow[6] = np.sum(marker_event_C5 == 1)
        else:
            print("Please input a number of cuts within 1 ~ 5.")
        print(jet_C1, bjet_C1, jet_C2, bjet_C2, jet_C3, bjet_C3, jet_C4, bjet_C4, jet_C5, bjet_C5)
        if num_of_cuts == 1:
            x = np.linspace(0,2,3)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5],[ "Origin", "C1"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75,"C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5".format(jet_C1, bjet_C1))
            plt.savefig(OUTPUT_FILE)
        elif num_of_cuts == 2:
            x = np.linspace(0,3,4)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5, 2.5],[ "Origin", "C1", "C2"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75,"C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5\nC2: {2} jets and {3} bjets passing pt>25 GeV and |eta| < 2.5\n".format(jet_C1, bjet_C1, jet_C2, bjet_C2))
            plt.savefig(OUTPUT_FILE)
        elif num_of_cuts == 3:
            x = np.linspace(0,4,5)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5, 2.5, 3.5],[ "Origin", "C1", "C2","C3"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75,"C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5\nC2: {2} jets and {3} bjets passing pt>25 GeV and |eta| < 2.5\nC3: {4} jets and {5} bjets passing pt>25 GeV and |eta| < 2.5\n".format(jet_C1, bjet_C1, jet_C2, bjet_C2, jet_C3, bjet_C3))
            plt.savefig(OUTPUT_FILE)
        elif num_of_cuts == 4:
            x = np.linspace(0,5,6)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5],[ "Origin", "C1", "C2", "C3", "C4"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75,"C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5\nC2: {2} jets and {3} bjets passing pt>25 GeV and |eta| < 2.5\nC3: {4} jets and {5} bjets passing pt>25 GeV and |eta| < 2.5\nC4: {6} jets with {7} bjet passing pt>25 GeV and |eta| < 2.5\n".format(jet_C1, bjet_C1, jet_C2, bjet_C2, jet_C3, bjet_C3, jet_C4, bjet_C4))
            plt.savefig(OUTPUT_FILE)
        elif num_of_cuts == 5:
            x = np.linspace(0,6,7)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5],[ "Origin", "C1", "C2", "C3", "C4", "C5"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75, "C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5\nC2: {2} jets and {3} bjets passing pt>25 GeV and |eta| < 2.5\nC3: {4} jets and {5} bjets passing pt>25 GeV and |eta| < 2.5\nC4: {6} jets with {7} bjet passing pt>25 GeV and |eta| < 2.5\nC5: {8} jet with {9} bjets passing pt>25 GeV and |eta| < 2.5".format(jet_C1, bjet_C1, jet_C2, bjet_C2, jet_C3, bjet_C3, jet_C4, bjet_C4, jet_C5, bjet_C5))
            plt.savefig(OUTPUT_FILE)
            
        else:
            print("Please input a number of cuts within 1 ~ 5.")
    
    elif int(SINGLE) != 1 and bool(SINGLE.isdigit) == True:
        _data = []
        particle = []
        jet = []
        files = os.listdir(INPUT_FILE)
        num_of_files = len(files)
        pbar = tqdm.tqdm(total=num_of_files)
        for a in files:
            try:
                print("Loading root file from {0}.".format(os.path.join(INPUT_FILE, a)))
                tmp = uproot.open(os.path.join(INPUT_FILE, a))['Delphes']
                _data.append(tmp)
                pbar.update(1) 
            except:
                print('Please check input file path.')
        
        for i in tqdm.trange(num_of_files):
            tmp_particle = particle_properties(_data[i])
            tmp_jet = jet_properties(_data[i])

            if i == 0 :
                _particle_event = tmp_particle.event
                _particle_pt = tmp_particle.pt
                _particle_eta = tmp_particle.eta
                _particle_phi = tmp_particle.phi
                _particle_pid = tmp_particle.pid
                _particle_M1 = tmp_particle.M1
                _particle_M2 = tmp_particle.M2
                _particle_D1 = tmp_particle.D1
                _particle_D2 = tmp_particle.D2
                _particle_status = tmp_particle.status
                _particle_rapidity = tmp_particle.rapidity
                _particle_mass = tmp_particle.mass
                _particle_charge = tmp_particle.charge

                _jet_event = tmp_jet.event
                _jet_pt = tmp_jet.pt
                _jet_eta = tmp_jet.eta
                _jet_phi = tmp_jet.phi
                _jet_btag = tmp_jet.btag
                _jet_area = tmp_jet.area
                _jet_mass = tmp_jet.mass
                _jet_charge = tmp_jet.charge
                
            else :
                _particle_event = np.concatenate((_particle_event, tmp_particle.event))
                _particle_pt = np.concatenate((_particle_pt, tmp_particle.pt))
                _particle_eta = np.concatenate((_particle_eta, tmp_particle.eta))
                _particle_phi = np.concatenate((_particle_phi, tmp_particle.phi))
                _particle_pid = np.concatenate((_particle_pid, tmp_particle.pid))
                _particle_M1 = np.concatenate((_particle_M1, tmp_particle.M1))
                _particle_M2 = np.concatenate((_particle_M2, tmp_particle.M2))
                _particle_D1 = np.concatenate((_particle_D1, tmp_particle.D1))
                _particle_D2 = np.concatenate((_particle_D2, tmp_particle.D2))
                _particle_status = np.concatenate((_particle_status, tmp_particle.status))
                _particle_rapidity = np.concatenate((_particle_rapidity, tmp_particle.rapidity))
                _particle_mass = np.concatenate((_particle_mass, tmp_particle.mass))
                _particle_charge = np.concatenate((_particle_charge, tmp_particle.charge))

                _jet_event = np.concatenate((_jet_event, tmp_jet.event))
                _jet_pt = np.concatenate((_jet_pt, tmp_jet.pt))
                _jet_eta = np.concatenate((_jet_eta, tmp_jet.eta))
                _jet_phi = np.concatenate((_jet_phi, tmp_jet.phi))
                _jet_btag = np.concatenate((_jet_btag, tmp_jet.btag))
                _jet_area = np.concatenate((_jet_area, tmp_jet.area))
                _jet_mass = np.concatenate((_jet_mass, tmp_jet.mass))
                _jet_charge = np.concatenate((_jet_charge, tmp_jet.charge))

        clustered_particle_data = properties_for_particle(_particle_event, _particle_pt, _particle_eta, _particle_phi, _particle_pid, _particle_M1, _particle_M2, _particle_D1, _particle_D2, _particle_status, _particle_rapidity, _particle_mass, _particle_charge)

        clustered_jet_data = properties_for_jet(_jet_event, _jet_pt, _jet_eta, _jet_phi, _jet_btag, _jet_area, _jet_mass, _jet_charge)
        
        if num_of_cuts == 1:
            marker_event_C1 = []
            marker_jet = []
            marker_bjet = []
            for i in range(len(clustered_particle_data.event)):
                marker_event_C1.append(0)
                marker_jet.append(np.zeros([len(clustered_jet_data.pt[i])]))
                marker_bjet.append(np.zeros([len(clustered_jet_data.pt[i])]))
            
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)

        elif num_of_cuts == 2:
            marker_event_C1 = []
            marker_event_C2 = []
            marker_jet = []
            marker_bjet = []
            for i in range(len(clustered_particle_data.event)):
                marker_event_C1.append(0)
                marker_event_C2.append(0)
                marker_jet.append(np.zeros([len(clustered_jet_data.pt[i])]))
                marker_bjet.append(np.zeros([len(clustered_jet_data.pt[i])]))
            
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_event_C2 = np.asanyarray(marker_event_C2, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)

        elif num_of_cuts == 3:
            marker_event_C1 = []
            marker_event_C2 = []
            marker_event_C3 = []
            marker_jet = []
            marker_bjet = []
            for i in range(len(clustered_particle_data.event)):
                marker_event_C1.append(0)
                marker_event_C2.append(0)
                marker_event_C3.append(0)
                marker_jet.append(np.zeros([len(clustered_jet_data.pt[i])]))
                marker_bjet.append(np.zeros([len(clustered_jet_data.pt[i])]))
            
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_event_C2 = np.asanyarray(marker_event_C2, dtype=object)
            marker_event_C3 = np.asanyarray(marker_event_C3, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)

        elif num_of_cuts == 4:
            marker_event_C1 = []
            marker_event_C2 = []
            marker_event_C3 = []
            marker_event_C4 = []
            marker_jet = []
            marker_bjet = []

            for i in range(len(clustered_particle_data.event)):
                marker_event_C1.append(0)
                marker_event_C2.append(0)
                marker_event_C3.append(0)
                marker_event_C4.append(0)
                marker_jet.append(np.zeros([len(clustered_jet_data.pt[i])]))
                marker_bjet.append(np.zeros([len(clustered_jet_data.pt[i])]))
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_event_C2 = np.asanyarray(marker_event_C2, dtype=object)
            marker_event_C3 = np.asanyarray(marker_event_C3, dtype=object)
            marker_event_C4 = np.asanyarray(marker_event_C4, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)
        elif num_of_cuts == 5:
            marker_event_C1 = []
            marker_event_C2 = []
            marker_event_C3 = []
            marker_event_C4 = []
            marker_event_C5 = []
            marker_jet = []
            marker_bjet = []
            for i in range(len(clustered_particle_data.event)):
                marker_event_C1.append(0)
                marker_event_C2.append(0)
                marker_event_C3.append(0)
                marker_event_C4.append(0)
                marker_event_C5.append(0)
                marker_jet.append(np.zeros([len(clustered_jet_data.pt[i])]))
                marker_bjet.append(np.zeros([len(clustered_jet_data.pt[i])]))
            marker_event_C1 = np.asanyarray(marker_event_C1, dtype=object)
            marker_event_C2 = np.asanyarray(marker_event_C2, dtype=object)
            marker_event_C3 = np.asanyarray(marker_event_C3, dtype=object)
            marker_event_C4 = np.asanyarray(marker_event_C4, dtype=object)
            marker_event_C5 = np.asanyarray(marker_event_C5, dtype=object)
            marker_jet = np.asanyarray(marker_jet, dtype=object)
            marker_bjet = np.asanyarray(marker_bjet, dtype=object)
        else:
            print("Please input a number of cuts within 1 ~ 5.")
        

        for i in tqdm.trange(len(clustered_particle_data.event)):
            for j in range(len(clustered_jet_data.pt[i])):
                if clustered_jet_data.btag[i][j] == 1 and clustered_jet_data.pt[i][j] > pt_cuts and np.abs(clustered_jet_data.eta[i][j]) < eta_cuts:
                    marker_bjet[i][j] = 1 
                else: pass 
            
                if clustered_jet_data.pt[i][j] > pt_cuts and np.abs(clustered_jet_data.eta[i][j]) <= eta_cuts:
                    marker_jet[i][j] = 1
                else: pass 

        for i in tqdm.trange(len(clustered_particle_data.event)):
            if num_of_cuts == 1:
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
            elif num_of_cuts == 2:
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C2 and np.sum(marker_bjet[i] == 1) >= bjet_C2:
                    marker_event_C2[i] = 1 
            elif num_of_cuts == 3:
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C2 and np.sum(marker_bjet[i] == 1) >= bjet_C2:
                    marker_event_C2[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C3 and np.sum(marker_bjet[i] == 1) >= bjet_C3:
                    marker_event_C3[i] = 1 
            elif num_of_cuts == 4: 
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C2 and np.sum(marker_bjet[i] == 1) >= bjet_C2:
                    marker_event_C2[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C3 and np.sum(marker_bjet[i] == 1) >= bjet_C3:
                    marker_event_C3[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C4 and np.sum(marker_bjet[i] == 1) >= bjet_C4:
                    marker_event_C4[i] = 1 
            elif num_of_cuts == 5:
                if np.sum(marker_jet[i] == 1) >= jet_C1 and np.sum(marker_bjet[i] == 1) >= bjet_C1:
                    marker_event_C1[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C2 and np.sum(marker_bjet[i] == 1) >= bjet_C2:
                    marker_event_C2[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C3 and np.sum(marker_bjet[i] == 1) >= bjet_C3:
                    marker_event_C3[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C4 and np.sum(marker_bjet[i] == 1) >= bjet_C4:
                    marker_event_C4[i] = 1 
                if np.sum(marker_jet[i] == 1) >= jet_C5 and np.sum(marker_bjet[i] == 1) >= bjet_C5:
                    marker_event_C5[i] = 1 

        if num_of_cuts == 1:
            cutflow = np.zeros(3)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
        elif num_of_cuts == 2:
            cutflow = np.zeros(4)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
            cutflow[3] = np.sum(marker_event_C2 == 1)
        elif num_of_cuts == 3:
            cutflow = np.zeros(5)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
            cutflow[3] = np.sum(marker_event_C2 == 1)
            cutflow[4] = np.sum(marker_event_C3 == 1)
        elif num_of_cuts == 4:
            cutflow = np.zeros(6)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
            cutflow[3] = np.sum(marker_event_C2 == 1)
            cutflow[4] = np.sum(marker_event_C3 == 1)
            cutflow[5] = np.sum(marker_event_C4 == 1)
        elif num_of_cuts == 5:
            cutflow = np.zeros(7)
            cutflow[0] = len(marker_event_C1)
            cutflow[1] = len(marker_event_C1)
            cutflow[2] = np.sum(marker_event_C1 == 1)
            cutflow[3] = np.sum(marker_event_C2 == 1)
            cutflow[4] = np.sum(marker_event_C3 == 1)
            cutflow[5] = np.sum(marker_event_C4 == 1)
            cutflow[6] = np.sum(marker_event_C5 == 1)
        else:
            print("Please input a number of cuts within 1 ~ 5.")

        if num_of_cuts == 1:
            x = np.linspace(0,2,3)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5],[ "Origin", "C1"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75,"C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5".format(jet_C1, bjet_C1))
            plt.savefig(OUTPUT_FILE)
        elif num_of_cuts == 2:
            x = np.linspace(0,3,4)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5, 2.5],[ "Origin", "C1", "C2"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75,"C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5\nC2: {2} jets and {3} bjets passing pt>25 GeV and |eta| < 2.5\n".format(jet_C1, bjet_C1, jet_C2, bjet_C2))
            plt.savefig(OUTPUT_FILE)
        elif num_of_cuts == 3:
            x = np.linspace(0,4,5)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5, 2.5, 3.5],[ "Origin", "C1", "C2","C3"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75,"C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5\nC2: {2} jets and {3} bjets passing pt>25 GeV and |eta| < 2.5\nC3: {4} jets and {5} bjets passing pt>25 GeV and |eta| < 2.5\n".format(jet_C1, bjet_C1, jet_C2, bjet_C2, jet_C3, bjet_C3))
            plt.savefig(OUTPUT_FILE)
        elif num_of_cuts == 4:
            x = np.linspace(0,5,6)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5],[ "Origin", "C1", "C2", "C3", "C4"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.75,"C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5\nC2: {2} jets and {3} bjets passing pt>25 GeV and |eta| < 2.5\nC3: {4} jets and {5} bjets passing pt>25 GeV and |eta| < 2.5\nC4: {6} jets with {7} bjet passing pt>25 GeV and |eta| < 2.5\n".format(jet_C1, bjet_C1, jet_C2, bjet_C2, jet_C3, bjet_C3, jet_C4, bjet_C4))
            plt.savefig(OUTPUT_FILE)
        elif num_of_cuts == 5:
            x = np.linspace(0,6,7)
            plt.figure(figsize=(8,6))
            plt.step(x, cutflow)
            plt.xlabel("Cuts")
            plt.ylabel("# of events")
            plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5],[ "Origin", "C1", "C2", "C3", "C4", "C5"])
            x0, xmax = plt.xlim()
            y0, ymax = plt.ylim()
            data_width = xmax - x0
            data_height = ymax - y0
            plt.text(data_width*0.6,data_height*0.5, "C1: {0} jets and {1} bjets passing pt>25 GeV and |eta| < 2.5\nC2: {2} jets and {3} bjets passing pt>25 GeV and |eta| < 2.5\nC3: {4} jets and {5} bjets passing pt>25 GeV and |eta| < 2.5\nC4: {6} jets with {7} bjet passing pt>25 GeV and |eta| < 2.5\nC5: {8} jet with {9} bjets passing pt>25 GeV and |eta| < 2.5".format(jet_C1, bjet_C1, jet_C2, bjet_C2, jet_C3, bjet_C3, jet_C4, bjet_C4, jet_C5, bjet_C5),horizontalalignment='center',verticalalignment='center')
            plt.savefig(OUTPUT_FILE)
            
        else:
            print("Please input a number of cuts within 1 ~ 5.")
        print("+------------------------------------------------------------------------------------------------------+")
        print("Figure has been saved. Path: {0}.".format(OUTPUT_FILE))
        print("+------------------------------------------------------------------------------------------------------+")
    else:
        print("Please inpurt a correct mode.\n1.-s 1\n2.-s > 1")
    
    

