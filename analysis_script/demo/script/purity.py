"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
import pandas as pd 
import numpy as np 
import h5py, sys, traceback, os, tqdm, time
from .utilize import pdgid, purity_classifier
import multiprocessing as mp



def purity(INPUT_FILE, TARGET_FILE, OUTPUT_FILE, MODEL, SINGLE):
    if SINGLE == 1:
        with np.load(TARGET_FILE, allow_pickle=True) as dalta_R:
            jet_parton_index_del = dalta_R['jet_parton_index'][:]
            jet_barcode_del = dalta_R['jet_barcode'][:]
            jet_pt_del = dalta_R['jet_pt'][:]
            jet_eta_del = dalta_R['jet_eta'][:]
            jet_phi_del = dalta_R['jet_phi'][:]
            jet_mass_del = dalta_R['jet_mass'][:]
            jet_btag_del = dalta_R['jet_btag'][:]
            parton_jet_index_del = dalta_R['parton_jet_index'][:]
            parton_pdgid_del = dalta_R['parton_pdgid'][:]
            parton_barcode_del = dalta_R['parton_barcode'][:]
            parton_pt_del = dalta_R['parton_pt'][:]
            parton_eta_del = dalta_R['parton_eta'][:]
            parton_phi_del = dalta_R['parton_phi'][:]
            parton_mass_del = dalta_R['parton_mass'][:]
            N_match_top_in_event_del = dalta_R['N_match_top_in_event'][:]
        with np.load(INPUT_FILE, allow_pickle=True) as f: 
            jet_parton_index_chi2 = f['jet_parton_index'][:]
            jet_barcode_chi2 = f['jet_barcode'][:]
            jet_pt_chi2 = f['jet_pt'][:]
            jet_eta_chi2 = f['jet_eta'][:]
            jet_phi_chi2 = f['jet_phi'][:]
            jet_mass_chi2 = f['jet_mass'][:]
            jet_btag_chi2 = f['jet_btag'][:]
            parton_jet_index_chi2 = f['parton_jet_index'][:]
            parton_pdgid_chi2 = f['parton_pdgid'][:]
            parton_barcode_chi2 = f['parton_barcode'][:]
            parton_pt_chi2 = f['parton_pt'][:]
            parton_eta_chi2 = f['parton_eta'][:]
            parton_phi_chi2 = f['parton_phi'][:]
            parton_mass_chi2 = f['parton_mass'][:]
            N_match_top_in_event_chi2 = f['N_match_top_in_event'][:]
    if SINGLE != 1:
        input_file_list = os.listdir(INPUT_FILE)
        target_file_list = os.listdif(TARGET_FILE)
        num_of_input_files = len(input_file_list)
        num_of_target_files = len(target_file_list)

        input_pbar = tqdm.tqdm(total=num_of_input_files)
        count = 0
        for a, b in zip(input_file_list, target_file_list):
            try:
                if count == 0:
                    print("Loading root file from {0}.".format(os.path.join(INPUT_FILE, a)))
                    with np.load(os.path.join(INPUT_FILE, a), allow_pickle=True) as f: 
                        jet_parton_index_chi2 = f['jet_parton_index'][:]
                        jet_barcode_chi2 = f['jet_barcode'][:]
                        jet_pt_chi2 = f['jet_pt'][:]
                        jet_eta_chi2 = f['jet_eta'][:]
                        jet_phi_chi2 = f['jet_phi'][:]
                        jet_mass_chi2 = f['jet_mass'][:]
                        jet_btag_chi2 = f['jet_btag'][:]
                        parton_jet_index_chi2 = f['parton_jet_index'][:]
                        parton_pdgid_chi2 = f['parton_pdgid'][:]
                        parton_barcode_chi2 = f['parton_barcode'][:]
                        parton_pt_chi2 = f['parton_pt'][:]
                        parton_eta_chi2 = f['parton_eta'][:]
                        parton_phi_chi2 = f['parton_phi'][:]
                        parton_mass_chi2 = f['parton_mass'][:]
                        N_match_top_in_event_chi2 = f['N_match_top_in_event'][:]

                    print("Loading root file from {0}.".format(os.path.join(TARGET_FILE, b)))
                    with np.load(os.path.join(TARGET_FILE, b), allow_pickle=True) as f: 
                        jet_parton_index_del = f['jet_parton_index'][:]
                        jet_barcode_del = f['jet_barcode'][:]
                        jet_pt_del = f['jet_pt'][:]
                        jet_eta_del = f['jet_eta'][:]
                        jet_phi_del = f['jet_phi'][:]
                        jet_mass_del = f['jet_mass'][:]
                        jet_btag_del = f['jet_btag'][:]
                        parton_jet_index_del = f['parton_jet_index'][:]
                        parton_pdgid_del = f['parton_pdgid'][:]
                        parton_barcode_del = f['parton_barcode'][:]
                        parton_pt_del = f['parton_pt'][:]
                        parton_eta_del = f['parton_eta'][:]
                        parton_phi_del = f['parton_phi'][:]
                        parton_mass_del = f['parton_mass'][:]
                        N_match_top_in_event_del = f['N_match_top_in_event'][:]
                    pbar.update(1) 
                else: 
                    print("Loading root file from {0}.".format(os.path.join(INPUT_FILE, a)))
                    with np.load(os.path.join(INPUT_FILE, a), allow_pickle=True) as f: 
                        jet_parton_index_chi2 = np.concatenate((jet_parton_index_chi2, f['jet_parton_index'][:]))
                        jet_barcode_chi2 = np.concatenate((jet_barcode_chi2, f['jet_barcode'][:]))
                        jet_pt_chi2 = np.concatenate((jet_pt_chi2, f['jet_pt'][:]))
                        jet_eta_chi2 = np.concatenate((jet_eta_chi2, f['jet_eta'][:]))
                        jet_phi_chi2 = np.concatenate((jet_phi_chi2, f['jet_phi'][:]))
                        jet_mass_chi2 = np.concatenate((jet_mass_chi2, f['jet_mass'][:]))
                        jet_btag_chi2 = np.concatenate((jet_btag_chi2, f['jet_btag'][:]))
                        parton_jet_index_chi2 = np.concatenate((parton_jet_index_chi2, f['parton_jet_index'][:]))
                        parton_pdgid_chi2 = np.concatenate((parton_pdgid_chi2, f['parton_pdgid'][:]))
                        parton_barcode_chi2 = np.concatenate((parton_barcode_chi2, f['parton_barcode'][:]))
                        parton_pt_chi2 = np.concatenate((parton_pt_chi2, f['parton_pt'][:]))
                        parton_eta_chi2 = np.concatenate((parton_eta_chi2, f['parton_eta'][:]))
                        parton_phi_chi2 = np.concatenate((parton_phi_chi2, f['parton_phi'][:]))
                        parton_mass_chi2 = np.concatenate((parton_mass_chi2, f['parton_mass'][:]))
                        N_match_top_in_event_chi2 = np.concatenate((N_match_top_in_event_chi2, f['N_match_top_in_event'][:]))

                    print("Loading root file from {0}.".format(os.path.join(TARGET_FILE, b)))
                    with np.load(os.path.join(TARGET_FILE, b), allow_pickle=True) as f: 
                        jet_parton_index_del = np.concatenate((jet_parton_index_del, f['jet_parton_index'][:]))
                        jet_barcode_del = np.concatenate((jet_barcode_del, f['jet_barcode'][:]))
                        jet_pt_del = np.concatenate((jet_pt_del, f['jet_pt'][:]))
                        jet_eta_del = np.concatenate((jet_eta_del, f['jet_eta'][:]))
                        jet_phi_del = np.concatenate((jet_phi_del, f['jet_phi'][:]))
                        jet_mass_del = np.concatenate((jet_mass_del, f['jet_mass'][:]))
                        jet_btag_del = np.concatenate((jet_btag_del, f['jet_btag'][:]))
                        parton_jet_index_del = np.concatenate((parton_jet_index_del, f['parton_jet_index'][:]))
                        parton_pdgid_del = np.concatenate((parton_pdgid_del, f['parton_pdgid'][:]))
                        parton_barcode_del = np.concatenate((parton_barcode_del, f['parton_barcode'][:]))
                        parton_pt_del = np.concatenate((parton_pt_del, f['parton_pt'][:]))
                        parton_eta_del = np.concatenate((parton_eta_del, f['parton_eta'][:]))
                        parton_phi_del = np.concatenate((parton_phi_del, f['parton_phi'][:]))
                        parton_mass_del = np.concatenate((parton_mass_del, f['parton_mass'][:]))
                        N_match_top_in_event_del = np.concatenate((N_match_top_in_event_del, f['N_match_top_in_event'][:]))
                    pbar.update(1) 
            except:
                print('Please check input file path.')
    print("+------------------------------------------------------------------------------------------------------+")
    print(f'Start computing purity, model: {MODEL}.')
    print("+------------------------------------------------------------------------------------------------------+")
    
    #Both correct ----> case 1
    #1 correct + 1 incorrect. ----> case 2
    #Neither correct. ----> case 3
    #1 or 2 Non-matched candidate exist. ----> case 4
    if MODEL == 'ttbar':
        correct_left = np.zeros(len(parton_jet_index_del))
        correct_right = np.zeros(len(parton_jet_index_del))
        case = np.zeros(len(parton_jet_index_del))

        for i in tqdm.trange(len(parton_barcode_del)):

            left_target = parton_jet_index_chi2[i][:3]
            right_target = parton_jet_index_chi2[i][3:]

            left_src = parton_jet_index_del[i][:3]    
            right_src = parton_jet_index_del[i][3:]  

            if np.sum(left_src == 'Nan') >0:
                correct_left[i] = 4       
            if np.sum(right_src == 'Nan') >0:
                correct_right[i] = 4

            if correct_left[i] != 4 and correct_right[i] != 4:
                correct_left[i], correct_right[i] = purity_classifier(parton_jet_index_del[i], parton_jet_index_chi2[i], "pair", MODEL)
            elif correct_left[i] != 4 and correct_right[i] == 4:
                correct_left[i] = purity_classifier(parton_jet_index_del[i], parton_jet_index_chi2[i], "left", MODEL)
            elif correct_left[i] == 4 and correct_right[i] != 4:
                correct_right[i] = purity_classifier(parton_jet_index_del[i], parton_jet_index_chi2[i], "right", MODEL)
            else:
                pass

            if correct_left[i] == 1 and correct_right[i] == 1:
                case[i] = 1
            elif correct_left[i] != 1 and correct_right[i] == 1:
                case[i] = 2
            elif correct_left[i] == 1 and correct_right[i] != 1:
                case[i] = 2
            elif  correct_left[i] != 1 and  correct_right[i] != 1 and correct_left[i] != 4 and  correct_right[i] != 4:
                case[i] = 3
            elif correct_left[i] == 4 or  correct_right[i] == 4:
                case[i] = 4
    else:
        print("Please select a available model. (Only ttbar model(fully hadronic decay) available currently.)")
    
    print("+------------------------------------------------------------------------------------------------------+")
    print('Finished to compute purity.')
    print("+------------------------------------------------------------------------------------------------------+")
    
    print("+------------------------------------------------------------------------------------------------------+")
    print("Writing event record to the npz file.")
    print("+------------------------------------------------------------------------------------------------------+")
    np.savez_compressed(OUTPUT_FILE, truth_match_result=parton_jet_index_del, chi2_reconstructresult=parton_jet_index_chi2, case=case)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Event record has been send to {0}.npz.".format(OUTPUT_FILE))
    print("+------------------------------------------------------------------------------------------------------+")