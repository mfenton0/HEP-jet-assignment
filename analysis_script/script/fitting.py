"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
import numpy as np
import pandas as pd
import h5py, os, tqdm, sys, uproot

def cal_two_jet_inv(jet1, jet2):
            part_1 = (jet1.e + jet2.e)**2
            part_2 = (jet1.px + jet2.px)**2
            part_3 = (jet1.py + jet2.py)**2
            part_4 = (jet1.pz + jet2.pz)**2
            return np.sqrt( part_1 - part_2 - part_3 - part_4 )

def cal_three_jet_inv(jet1, jet2, jet3):
    part_1 = (jet1.e + jet2.e + jet3.e)**2
    part_2 = (jet1.px + jet2.px + jet3.px)**2
    part_3 = (jet1.py + jet2.py + jet3.py)**2
    part_4 = (jet1.pz + jet2.pz + jet3.pz)**2
    return np.sqrt( part_1 - part_2 - part_3 - part_4 )

def fitting(INPUT_FILE, OUTPUT_FILE, MODEL, SINGLE):

    if SINGLE:
        try:
            print("Loading npz file from {0}.".format(INPUT_FILE))
            
            with np.load(INPUT_FILE, allow_pickle=True) as f:
                jet_parton_index = f['jet_parton_index'][:]
                parton_jet_index = f['parton_jet_index'][:]
                N_match_top_in_event = f['N_match_top_in_event'][:]
                if MODEL == "ttH":
                    N_match_higgs_in_event = f['N_match_higgs_in_event'][:]
                jet_pt = f['jet_pt'][:]
                jet_eta = f['jet_eta'][:]
                jet_phi = f['jet_phi'][:]
                jet_mass = f['jet_mass'][:]
        except:
            print('Please check input file path.')

    else :
        files = os.listdir(INPUT_FILE)
        num_of_files = len(files)
        pbar = tqdm.tqdm(total=num_of_files)
        for i in range(len(files)):
            if i == 0:
                try:
                    print("Loading hdf5 file from {0}.".format(files[i]))
                    with np.load(INPUT_FILE, allow_pickle=True) as f:
                        jet_parton_index = f['jet_parton_index'][:]
                        parton_jet_index = f['parton_jet_index'][:]
                        N_match_top_in_event = f['N_match_top_in_event'][:]
                        if MODEL == "ttH":
                            N_match_higgs_in_event = f['N_match_higgs_in_event'][:]
                        jet_pt = f['jet_pt'][:]
                        jet_eta = f['jet_eta'][:]
                        jet_phi = f['jet_phi'][:]
                        jet_mass = f['jet_mass'][:]
                    pbar.update(1) 
                except:
                    print('Please check input file path.')
            else : 
                try:
                    print("Loading root file from {0}.".format(files[i]))
                    with np.load(INPUT_FILE, allow_pickle=True) as f:
                        jet_parton_index = np.concatenate((jet_parton_index, f['jet_parton_index'][:]))
                        parton_jet_index = np.concatenate((parton_jet_index, f['parton_jet_index'][:]))
                        N_match_top_in_event = np.concatenate((N_match_top_in_event, f['N_match_top_in_event'][:]))
                        if MODEL == "ttH":
                            N_match_higgs_in_event = np.concatenate((N_match_higgs_in_event, f['N_match_higgs_in_event'][:]))
                        jet_pt = np.concatenate((jet_pt, f['jet_pt'][:]))
                        jet_eta = np.concatenate((jet_eta, f['jet_eta'][:]))
                        jet_phi = np.concatenate((jet_phi, f['jet_phi'][:]))
                        jet_mass =np.concatenate((jet_mass, f['jet_mass'][:]))
                    pbar.update(1) 
                except:
                    print('Please check input file path.')

    class candidate_property():
        def __init__(self, index_1, index_2):
            self.pt = jet_pt[index_1][index_2]
            self.eta = jet_eta[index_1][index_2]
            self.phi = jet_phi[index_1][index_2]
            self.mass = jet_mass[index_1][index_2]
            self.px = self.pt*np.cos(self.phi) 
            self.py = self.pt*np.sin(self.phi)
            self.pz = self.pt*np.sinh(self.eta)
            self.e = np.sqrt( (self.px**2 + self.py**2 + self.pz**2) + self.mass**2 )

    W_inv = []
    W_minus_inv = []
    top_inv = []
    top_bar_inv =  []
    

    for i in range(len(jet_parton_index)):
    # for i in range(10):
        if N_match_top_in_event[i] == 2:
            # print(jet_parton_index[i], N_match_top_in_event[i], N_match_higgs_in_event[i])
            for j in range(len(jet_parton_index[i])):
                if jet_parton_index[i][j] == 0:
                    # print(0)
                    _bjet_1 = candidate_property(i, j)
                elif jet_parton_index[i][j] == 1:
                    # print(1)
                    _jet_1 = candidate_property(i, j)
                elif jet_parton_index[i][j] == 2:
                    # print(2)
                    _jet_2 = candidate_property(i, j)
                elif jet_parton_index[i][j] == 3:
                    # print(3)
                    _bjet_2 = candidate_property(i, j)
                elif jet_parton_index[i][j] == 4:
                    # print(4)
                    _jet_3 = candidate_property(i, j)
                elif jet_parton_index[i][j] == 5:
                    # print(5)
                    _jet_4 = candidate_property(i, j)
                else : pass
        
            W_inv.append(cal_two_jet_inv(_jet_1, _jet_2))
            W_minus_inv.append(cal_two_jet_inv(_jet_3, _jet_4))
            top_inv.append(cal_three_jet_inv(_bjet_1, _jet_1, _jet_2))
            top_bar_inv.append(cal_three_jet_inv(_bjet_2, _jet_3, _jet_4))

    mean_W, sigma_W = np.mean(W_inv), np.std(W_inv)
    mean_W_minus, sigma_W_minus = np.mean(W_minus_inv), np.std(W_minus_inv)
    mean_top, sigma_top = np.mean(top_inv), np.std(top_inv)
    mean_top_bar, sigma_top_bar = np.mean(top_bar_inv), np.std(top_bar_inv)

    if MODEL == "ttH":
        higgs_inv = []
        for i in range(len(jet_parton_index)):
            if N_match_top_in_event[i] == 2 and N_match_higgs_in_event[i] == 1:
                for j in range(len(jet_parton_index[i])):
                    if jet_parton_index[i][j] == 6:
                        # print(6)
                        _bjet_3 = candidate_property(i, j)
                    elif jet_parton_index[i][j] == 7:
                        # print(7)
                        _bjet_4 = candidate_property(i, j)
                    else: pass
                    
                higgs_inv.append(cal_two_jet_inv(_bjet_3, _bjet_4))
        
        mean_h, sigma_h = np.mean(higgs_inv), np.std(higgs_inv)

    with open(OUTPUT_FILE, "w") as f:
        if MODEL == "ttbar":
            f.write("Mean value of W+ invariant mass is: {0}, sigma is: {1}\n".format(mean_W, sigma_W))
            f.write("Mean value of W- invariant mass is: {0}, sigma is: {1}\n".format(mean_W_minus, sigma_W_minus))
            f.write("Mean value of t invariant mass is: {0}, sigma is: {1}\n".format(mean_top, sigma_top))
            f.write("Mean value of t_bar invariant mass is: {0}, sigma is: {1}\n".format(mean_top_bar, sigma_top_bar))
        elif MODEL == "ttH":
            f.write("Mean value of W+ invariant mass is: {0}, sigma is: {1}\n".format(mean_W, sigma_W))
            f.write("Mean value of W- invariant mass is: {0}, sigma is: {1}\n".format(mean_W_minus, sigma_W_minus))
            f.write("Mean value of t invariant mass is: {0}, sigma is: {1}\n".format(mean_top, sigma_top))
            f.write("Mean value of t_bar invariant mass is: {0}, sigma is: {1}\n".format(mean_top_bar, sigma_top_bar))
            f.write("Mean value of Higgs boson mass is: {0}, sigma is: {1}".format(mean_h, sigma_h))
        elif MODEL == "four_top":
            print("Work in progress.")
        else: 
            print("Please select a correct model.")
        
