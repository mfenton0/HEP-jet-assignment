"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
#Import packages
import uproot
import pandas as pd 
import numpy as np 
from script.particle import particle_properties  #import particle properties helper function from particle_properties.py
from script.jet import jet_properties  #import jet properties helper function from jet_properties.py
from script.MissingET import Missing_ET_properties
from script.electron import electron_properties
from script.muon import muon_properties
import h5py, sys, traceback, os, tqdm, time
from script.utilize import delta_R, deltaPhi, pdgid, event_selection, quark_finder, deltaPhi, particle_tracing, chi_square_minimizer
import multiprocessing as mp

def chi2_from_npz(INPUT_FILE, OUTPUT_FILE, MODEL, SINGLE, PROCESS,EXTRA):
    if MODEL == "ttbar":
        """
        Barcode system
        t t~ W+ W- b b~ 
        0 0  0  0  0 0
        daughter of top and W+: 101000 ----> 40
        daughter of top and b: 101000 ----> 34
        daughter of anti top and W-: 100100 ----> 20
        daughter of anti top and b~: 100001 ----> 17
        """
        barcode = np.array([34, 40, 40, 17, 20, 20])
        NUM_OF_PARTON = 6
        NUM_OF_DAUGHTER = 6
    elif MODEL == "ttbar_lep_left":
        """
        Barcode system
        t t~ W+ W- b b~ 
        0 0  0  0  0 0
        daughter of top and W+: 101000 ----> 40
        daughter of top and b: 101000 ----> 34
        daughter of anti top and W-: 100100 ----> 20
        daughter of anti top and b~: 100001 ----> 17
        """
        barcode = np.array([34, 40, 40, 17, 20, 20])
        NUM_OF_PARTON = 4
        NUM_OF_DAUGHTER = 6
    elif MODEL == "ttbar_lep_right":
        """
        Barcode system
        t t~ W+ W- b b~ 
        0 0  0  0  0 0
        daughter of top and W+: 101000 ----> 40
        daughter of top and b: 101000 ----> 34
        daughter of anti top and W-: 100100 ----> 20
        daughter of anti top and b~: 100001 ----> 17
        """
        barcode = np.array([34, 40, 40, 17, 20, 20])
        NUM_OF_PARTON = 4
        NUM_OF_DAUGHTER = 6
    elif MODEL == "ttH":
        """
        Barcode system
        t t~ W+ W- b b~ H
        0 0  0  0  0 0  0
        daughter of t and b = 1000100  ----> 68
        daughter of t and W+ = 1010000 ----> 80
        daughter of t~ and W- = 0101000 ----> 34
        daughter of t~ and b~ = 0100010 ----> 40
        daughter of H = 0000001 ----> 1
        """
        barcode = np.array([68, 80, 80, 34, 40, 40, 1, 1])
        NUM_OF_PARTON = 8
        NUM_OF_DAUGHTER = 8
    elif MODEL == "four_top":
        """
        Barcode system
        t1 t2 t1~ t2~ W+1 W-1 W+2 W-2 b1 b2 b1~ b2~             describe          barcode   sequence
        0  0   0   0   0   0   0   0  0  0   0   0

        1  0   0   0   1   0   0   0  0  0   0   0  <--- daughter of t1 and W+1   2176         2,3
        1  0   0   0   0   0   0   0  1  0   0   0  <--- daughter of t1 and b1    2056          1
        0  0   1   0   0   1   0   0  0  0   0   0  <--- daughter of t1~ and W-1  576          5,6
        0  0   1   0   0   0   0   0  0  1   0   0  <--- daughter of t1~ and b1~  516           4

        0  1   0   0   0   0   1   0  0  0   0   0  <--- daughter of t2 and W+2   1056         7,8
        0  1   0   0   0   0   0   0  0  1   0   0  <--- daughter of t2 and b2    1028          9
        0  0   0   1   0   0   0   1  0  0   0   0  <--- daughter of t2~ and W-2  272          11.12
        0  0   0   1   0   0   0   0  0  0   0   1  <--- daughter of t2~ and b2~  257           10

        """
        barcode = np.array([2056, 2176, 2176, 516, 576, 576, 1028, 1056, 1056, 257, 272, 272])
        NUM_OF_PARTON = 12
        NUM_OF_DAUGHTER = 12
    else:
        print("Please select a correct model.")
        
    if int(SINGLE) == 1:
        with np.load(INPUT_FILE, allow_pickle=True) as file:
            jet_parton_index=file['jet_parton_index'][:]
            jet_barcode=file['jet_barcode'][:]
            jet_pt=file['jet_pt'][:]
            jet_eta=file['jet_eta'][:]
            jet_phi=file['jet_phi'][:]
            jet_mass=file['jet_mass'][:]
            jet_btag=file['jet_btag'][:]
            jet_num_of_jets=file['jet_num_of_jets'][:]
            parton_jet_index=file['parton_jet_index'][:]
            parton_pdgid=file['parton_pdgid'][:]
            parton_barcode=file['parton_barcode'][:]
            parton_pt=file['parton_pt'][:]
            parton_eta=file['parton_eta'][:]
            parton_phi=file['parton_phi'][:]
            parton_mass=file['parton_mass'][:]
            N_match_top_in_event=file['N_match_top_in_event'][:]
    elif int(SINGLE) != 1 and bool(SINGLE.isdigit) == True:
        files = os.listdir(INPUT_FILE)
        num_of_files = len(files)
        pbar = tqdm.tqdm(total=num_of_files)
        for i in range(len(files)):
            try:
                if i == 0:
                    print("Loading root file from {0}.".format(os.path.join(INPUT_FILE, files[i])))
                    with np.load(os.path.join(INPUT_FILE, files[i]), allow_pickle=True) as file:
                        jet_parton_index=file['jet_parton_index'][:]
                        jet_barcode=file['jet_barcode'][:]
                        jet_pt=file['jet_pt'][:]
                        jet_eta=file['jet_eta'][:]
                        jet_phi=file['jet_phi'][:]
                        jet_mass=file['jet_mass'][:]
                        jet_btag=file['jet_btag'][:]
                        jet_num_of_jets=file['jet_num_of_jets'][:]
                        parton_jet_index=file['parton_jet_index'][:]
                        parton_pdgid=file['parton_pdgid'][:]
                        parton_barcode=file['parton_barcode'][:]
                        parton_pt=file['parton_pt'][:]
                        parton_eta=file['parton_eta'][:]
                        parton_phi=file['parton_phi'][:]
                        parton_mass=file['parton_mass'][:]
                        N_match_top_in_event=file['N_match_top_in_event'][:]
                    pbar.update(1) 
                else:
                    print("Loading root file from {0}.".format(os.path.join(INPUT_FILE, files[i])))
                    with np.load(os.path.join(INPUT_FILE, files[i]), allow_pickle=True) as file:
                        jet_parton_index=np.concatenate((jet_parton_index, file['jet_parton_index'][:]))
                        jet_barcode=np.concatenate((jet_barcode,file['jet_barcode'][:]))
                        jet_pt=np.concatenate((jet_pt,file['jet_pt'][:]))
                        jet_eta=np.concatenate((jet_eta,file['jet_eta'][:]))
                        jet_phi=np.concatenate((jet_phi,file['jet_phi'][:]))
                        jet_mass=np.concatenate((jet_mass,file['jet_mass'][:]))
                        jet_btag=np.concatenate((jet_btag,file['jet_btag'][:]))
                        jet_num_of_jets=np.concatenate((jet_num_of_jets,file['jet_num_of_jets'][:]))
                        parton_jet_index=np.concatenate((parton_jet_index,file['parton_jet_index'][:]))
                        parton_pdgid=np.concatenate((parton_pdgid,file['parton_pdgid'][:]))
                        parton_barcode=np.concatenate((parton_barcode,file['parton_barcode'][:]))
                        parton_pt=np.concatenate((parton_pt,file['parton_pt'][:]))
                        parton_eta=np.concatenate((parton_eta,file['parton_eta'][:]))
                        parton_phi=np.concatenate((parton_phi,file['parton_phi'][:]))
                        parton_mass=np.concatenate((parton_mass,file['parton_mass'][:]))
                        N_match_top_in_event=np.concatenate((N_match_top_in_event,file['N_match_top_in_event'][:]))
                    
                    pbar.update(1)
            except:
                print('Please check input file path.')
    else: 
        print("Please input a vaild 'SINGLE' argument.")
        
    print("+------------------------------------------------------------------------------------------------------+")
    print("Starting chi-square matching.")
    print("+------------------------------------------------------------------------------------------------------+")
    start = time.time()

    _src_chi2 = [list([jet_pt[i], jet_eta[i], jet_phi[i], jet_btag[i], jet_mass[i], MODEL, EXTRA, jet_num_of_jets[i]]) for i in range(len(jet_pt))]
    
    print("Using {0} process for accelerating speed.".format(PROCESS))
    with mp.Pool(PROCESS) as p:
        _result_chi2 = p.starmap(chi_square_minimizer, _src_chi2)
        p.close()
        p.join()
    _result_chi2 = np.array(_result_chi2)

    min_chi2_value = _result_chi2[:, 0]
    parton_jet_index_chi2 = np.array([ x for x in _result_chi2[:, 1]])
    jet_parton_index_chi2 = np.array([ x for x in _result_chi2[:, 2]])
    smallest_10_chi2_candidate = np.array([ x for x in _result_chi2[:, 3]])
    smallest_10_chi2_value = np.array([ x for x in _result_chi2[:, 4]])
    if (min_chi2_value == -1).sum() != 0: print("There exist some events failed to compute chi-square reconstuction.")
        
    print("+------------------------------------------------------------------------------------------------------+")
    print("Chi-square matching finished. Cost: {0:.1f} s".format(time.time() - start))
    print("+------------------------------------------------------------------------------------------------------+")
    
    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording barcode information.")
    print("+------------------------------------------------------------------------------------------------------+")
    
    jet_barcode_chi2 = jet_parton_index_chi2.copy()

    for i in range(len(barcode)):
        jet_barcode_chi2 = np.where(jet_barcode_chi2 == i, barcode[i], jet_barcode_chi2) 

    print("+------------------------------------------------------------------------------------------------------+")
    print("Barcode information has beed record.")
    print("+------------------------------------------------------------------------------------------------------+")
    
    if MODEL == 'ttH':
        target = [i for i in range(NUM_OF_PARTON)]
        N_match_top_in_event_chi2 = np.zeros([len(jet_pt)])
        for i in tqdm.trange(len(jet_parton_index_chi2)):
            
            intersetion = set(target).intersection(jet_parton_index_chi2[i])
            N_match_top_in_event_chi2[i] = len(intersetion) // 3

        N_match_higgs_in_event_chi2 = np.zeros([len(jet_pt)])
        for i in range(len(jet_parton_index_chi2)):
            if 7 in jet_parton_index_chi2[i]:
                if 6 in jet_parton_index_chi2[i]:
                    N_match_higgs_in_event_chi2[i] = 1

        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2)  ))
        print("+------------------------------------------------------------------------------------------------------+")
    elif MODEL == 'ttbar':
        target = [i for i in range(NUM_OF_PARTON)]
        N_match_top_in_event_chi2 = np.zeros([len(jet_pt)])
        for i in tqdm.trange(len(jet_parton_index_chi2)):
            
            intersetion = set(target).intersection(jet_parton_index_chi2[i])
            N_match_top_in_event_chi2[i] = len(intersetion) // 3
        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2)  ))
        print("+------------------------------------------------------------------------------------------------------+")
    elif MODEL == 'four_top':
        target = [i for i in range(NUM_OF_PARTON)]
        N_match_top_in_event_chi2 = np.zeros([len(jet_pt)])
        for i in tqdm.trange(len(jet_parton_index_chi2)):
            
            intersetion = set(target).intersection(jet_parton_index_chi2[i])
            count_inter = 0
            if intersetion.intersection(set([0, 1, 2])) == {0,1,2}:
                count_inter += 1
            if intersetion.intersection(set([3, 4, 5])) == {3,4,5}:
                count_inter += 1
            if intersetion.intersection(set([6, 7, 8])) == {6,7,8}:
                count_inter += 1
            if intersetion.intersection(set([9, 10, 11])) == {9,10,11}:
                count_inter += 1

            N_match_top_in_event_chi2[i] = count_inter
        print("+------------------------------------------------------------------------------------------------------+")
        print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.\nFound {2} events with 3 ttbar candidate exist.\nFound {3} events with 4 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2), np.sum(N_match_top_in_event == 3), np.sum(N_match_top_in_event == 4)  ))
        print("+------------------------------------------------------------------------------------------------------+")
#     elif MODEL == 'ttbar_lep_left' or MODEL == 'ttbar_lep_right':
#         N_match_top_in_event_chi2 = np.zeros([len(jet_pt)])
#         target = [i for i in range(NUM_OF_PARTON)]
#         for i in tqdm.trange(len(jet_parton_index_chi2)):
#             intersetion = set(target).intersection(jet_parton_index_chi2[i])
#             if MET_delta_R_result[i] == 1 and lepton_delta_R_result[i] == 1:              
#                 if len(intersetion) == 4:
#                     N_match_top_in_event[i] = 2
#                 elif 3 in intersetion:
#                     N_match_top_in_event[i] = 1
#                 elif intersetion.intersection(set([0, 1, 2])) == {0,1,2} and (3 in intersetion) == False:
#                     N_match_top_in_event[i] = 1
#             else: 
#                 if intersetion.intersection(set([0, 1, 2])) == {0, 1, 2}:
#                     N_match_top_in_event[i] = 1
#                 else: pass
#         print("+------------------------------------------------------------------------------------------------------+")
#         print("Jet-parton matching section complete.\nFound {0} events with 1 ttbar candidate exist.\nFound {1} events with 2 ttbar candidate exist.".format( np.sum(N_match_top_in_event == 1), np.sum(N_match_top_in_event == 2) ))
        print("+------------------------------------------------------------------------------------------------------+")
    else : pass
    
    # Save selected events
    print("+------------------------------------------------------------------------------------------------------+")
    print("Writing event record to the npz file.")
    print("+------------------------------------------------------------------------------------------------------+")

    if MODEL == 'ttbar' or MODEL == 'four_top':
        np.savez_compressed(OUTPUT_FILE, 
                            jet_parton_index=jet_parton_index,
                            jet_parton_index_chi2=jet_parton_index_chi2,
                            jet_barcode=jet_barcode,
                            jet_barcode_chi2=jet_barcode_chi2,
                            jet_pt=jet_pt,
                            jet_eta=jet_eta,
                            jet_phi=jet_phi,
                            jet_mass=jet_mass,
                            jet_btag=jet_btag,
                            parton_jet_index=parton_jet_index,
                            parton_jet_index_chi2=parton_jet_index_chi2,
                            parton_pdgid=parton_pdgid,
                            parton_barcode=parton_barcode,
                            parton_pt=parton_pt,
                            parton_eta=parton_eta,
                            parton_phi=parton_phi,
                            parton_mass=parton_mass,
                            N_match_top_in_event=N_match_top_in_event,
                            N_match_top_in_event_chi2=N_match_top_in_event_chi2,
                            min_chi2_value=min_chi2_value, 
                            smallest_10_chi2_candidate=smallest_10_chi2_candidate,
                            smallest_10_chi2_value=smallest_10_chi2_value)
                           
#     elif MODEL == 'ttH':
#         np.savez_compressed(OUTPUT_FILE, 
#                             jet_parton_index=jet_parton_index,
#                             jet_barcode=jet_barcode,
#                             jet_pt=jet_pt,
#                             jet_eta=jet_eta,
#                             jet_phi=jet_phi,
#                             jet_mass=jet_mass,
#                             jet_btag=jet_btag,
#                             parton_jet_index=parton_jet_index,
#                             parton_pdgid=parton_pdgid,
#                             parton_barcode=parton_barcode,
#                             parton_pt=parton_pt,
#                             parton_eta=parton_eta,
#                             parton_phi=parton_phi,
#                             parton_mass=parton_mass,
#                             N_match_top_in_event=N_match_top_in_event,
#                             N_match_higgs_in_event=N_match_higgs_in_event)
#     elif MODEL == 'ttbar_lep_left' or MODEL == 'ttbar_lep_right':
#         np.savez_compressed(OUTPUT_FILE, 
#                             jet_parton_index=jet_parton_index,
#                             jet_barcode=jet_barcode,
#                             jet_pt=jet_pt,
#                             jet_eta=jet_eta,
#                             jet_phi=jet_phi,
#                             jet_mass=jet_mass,
#                             jet_btag=jet_btag,
#                             parton_jet_index=parton_jet_index,
#                             parton_pdgid=parton_pdgid,
#                             parton_barcode=parton_barcode,
#                             parton_pt=parton_pt,
#                             parton_eta=parton_eta,
#                             parton_phi=parton_phi,
#                             parton_mass=parton_mass,
#                             N_match_top_in_event=N_match_top_in_event,
#                             lepton_pt=lepton_pt,
#                             lepton_eta=lepton_eta,
#                             lepton_phi=lepton_phi,
#                             lepton_pdgid=lepton_pdgid,
#                             MET=MET,
#                             MET_ETA=MET_ETA,
#                             MET_PHI=MET_PHI,
#                             simulation_neutrino_pt=simulation_neutrino_pt,
#                             simulation_neutrino_eta=simulation_neutrino_eta,
#                             simulation_neutrino_phi=simulation_neutrino_phi,
#                             simulation_neutrino_pdgid=simulation_neutrino_pdgid,
#                             simulation_neutrino_barcode=simulation_neutrino_barcode,
#                             simulation_neutrino_mass=simulation_neutrino_mass,
#                             simulation_lepton_pt=simulation_lepton_pt,
#                             simulation_lepton_eta=simulation_lepton_eta,
#                             simulation_lepton_phi=simulation_lepton_phi,
#                             simulation_lepton_pdgid=simulation_lepton_pdgid,
#                             simulation_lepton_mass=simulation_lepton_mass,
#                             simulation_lepton_barcode=simulation_lepton_barcode)

    print("+------------------------------------------------------------------------------------------------------+")
    print("Event record has been send to {0}.npz.".format(OUTPUT_FILE))
    print("+------------------------------------------------------------------------------------------------------+")


    
