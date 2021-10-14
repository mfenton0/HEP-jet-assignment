"""
Author: David Ho^1, Hideki Okawa^2
Institute1: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Institute2: Fudan University, Shanghai, China
Mail: davidho@gapp.nthu.edu.tw, hideki.okawa@cern.ch
"""
#Import packages
import uproot
import pandas as pd 
import numpy as np 
import numba as nb
import h5py, sys, traceback, os, tqdm, time
from script.utilize import IO_module, chi_square_minimizer, helper, pdgid, deltaR_matching, process_methods, delta_R, deltaPhi
import multiprocessing as mp
from collections import OrderedDict
import faulthandler
faulthandler.enable()
def parse(INPUT_FILE, OUTPUT_FILE, MODEL, PROCESS, GENERATOR, SINGLE=True, COMPUTE_CHI2=False, **kargs):

    PID = pdgid()
    require_lepton = ["ttbar_lep", "ttbar_lep_left", "ttbar_lep_right"]
    # Setting `STATUS_CODE` for different shower generator.
    if GENERATOR == 'py8':
        STATUS_CODE = 22
    elif GENERATOR == 'herwig7':
        STATUS_CODE = 11
    else: 
        print("Please select a correct shower generator. 1. py8, 2. herwig7.")

    MAX_NUM_OF_JETS = 30

    # Setting barcode, `NUM_OF_PARTON`, and `NUM_OF_DAUGHTER` for different model
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
    elif MODEL == "ttbar_lep":
        """
        Barcode system
        t t~ W+ W- b b~
        0 0  0  0  0 0
        daughter of top/anti-top and leptonic W: 101000 ----> 40
        b on the leptonic side: 101000 ----> 34
        daughter of top/anti-top and hadronic W: 100100 ----> 20
        b on the hadronic side: 100001 ----> 17
        """
        barcode = np.array([34, 40, 40, 17, 20, 20])
        NUM_OF_PARTON = 4
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
    elif MODEL == "ZH":
        """
        Barcode system
        Z H W+ W- b b~ 
        0 0  0  0  0 0
        daughter of higgs and W+: 011000 ----> 24
        daughter of higgs and W-: 010100 ----> 20
        daughter of Z and b: 100010 ----> 34
        daughter of Z and b~: 100001 ----> 33
        """
        barcode = np.array([24, 24, 20, 20, 34, 33])
        NUM_OF_PARTON = 6
        NUM_OF_DAUGHTER = 6
    else:
        print("Please select a correct model.")

    print("+------------------------------------------------------------------------------------------------------+")
    print("Start loading dataset.")
    print("+------------------------------------------------------------------------------------------------------+")
    read_dataset = IO_module(INPUT_FILE, MODEL, MULTI = False)
    dataset = read_dataset.read_ROOT()
    print("+------------------------------------------------------------------------------------------------------+")
    print("Dataset loaded.")
    print("+------------------------------------------------------------------------------------------------------+")

    print("+------------------------------------------------------------------------------------------------------+")
    print("Start event selection.")
    print("+------------------------------------------------------------------------------------------------------+")
    if MODEL in require_lepton:
        marker_event = process_methods.event_selection(MODEL, 
                                        pt=dataset["jet"]["pt"], 
                                        eta=dataset["jet"]["eta"], 
                                        phi=dataset["jet"]["phi"],
                                        btag=dataset["jet"]["btag"],
                                        electron_pt=dataset["electron"]["pt"],
                                        electron_eta=dataset["electron"]["eta"],
                                        electron_phi=dataset["electron"]["phi"],
                                        muon_pt=dataset["muon"]["pt"],
                                        muon_eta=dataset["muon"]["eta"],
                                        muon_phi=dataset["muon"]["phi"],
                                        )
    else:
        marker_event = process_methods.event_selection(MODEL, pt=dataset["jet"]["pt"], eta=dataset["jet"]["eta"], btag=dataset["jet"]["btag"])
    passed = np.where(marker_event == 1)[0]
    print("+------------------------------------------------------------------------------------------------------+")
    print("Jet selection done. {0} events has been selected.".format(len(passed)))
    print("+------------------------------------------------------------------------------------------------------+")

    print("+------------------------------------------------------------------------------------------------------+")
    print("Starting parton tracing and looking for its daughter.")
    print("+------------------------------------------------------------------------------------------------------+")
    #Particle tracing and daughter finding section
    if int(PROCESS) == 1:

        if MODEL == 'ttbar_lep':
            daughter_t1 = [process_methods.lephad_finder(helper.to_dataframe(dataset["particle"], i), PID.top, 1, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of leptonic top quark.')]
            daughter_t2 = [process_methods.lephad_finder(helper.to_dataframe(dataset["particle"], i), PID.top, 2, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of hadronic top quark.')] # we don't care the sign of PID for lephad_finder
            daughter_t1_W = [process_methods.lephad_finder(helper.to_dataframe(dataset["particle"], passed[i]), PID.w_plus, 1, MODEL) for i in tqdm.trange(len(passed), desc='Finding daughters of leptonic W boson.')]
            daughter_t2_W = [process_methods.lephad_finder(helper.to_dataframe(dataset["particle"], passed[i]), PID.w_plus, 2, MODEL) for i in tqdm.trange(len(passed), desc='Finding daughters of hadronic W boson.')] # we don't care the sign of PID for lephad_finder
            daughter_t1_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in daughter_t1_W])
            daughter_t2_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in daughter_t2_W])

            daughter_t1_W_1 = daughter_t1_W_idx[:,0]
            daughter_t1_W_2 = daughter_t1_W_idx[:,1]
            daughter_t1_b = np.array([a["daughter_2_idx"] for a in daughter_t1])

            daughter_t2_W_1 = daughter_t2_W_idx[:,0]
            daughter_t2_W_2 = daughter_t2_W_idx[:,1]
            daughter_t2_b = np.array([a["daughter_1_idx"] for a in daughter_t2])

 
        elif MODEL == 'ttbar' or MODEL == 'ttbar_lep_left' or MODEL =='ttbar_lep_right':
            daughter_t1 = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], i), PID.top, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of top quark.')]
            daughter_t2 = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], i), PID.anti_top, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of anti-top quark.')]
            daughter_t1_W = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], passed[i]), PID.w_plus, MODEL) for i in tqdm.trange(len(passed), desc='Finding daughters of W+ boson.')]
            daughter_t2_W = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], passed[i]), PID.w_minus, MODEL) for i in tqdm.trange(len(passed), desc='Finding daughters of W- boson.')]
            daughter_t1_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in daughter_t1_W])
            daughter_t2_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in daughter_t2_W])

            daughter_t1_W_1 = daughter_t1_W_idx[:,0]
            daughter_t1_W_2 = daughter_t1_W_idx[:,1]
            daughter_t1_b = np.array([a["daughter_2_idx"] for a in daughter_t1])

            daughter_t2_W_1 = daughter_t2_W_idx[:,0]
            daughter_t2_W_2 = daughter_t2_W_idx[:,1]
            daughter_t2_b = np.array([a["daughter_1_idx"] for a in daughter_t2])

        elif MODEL == 'ttH':
            daughter_t1 = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], i), PID.top, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of top quark.')]
            daughter_t2 = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], i), PID.anti_top, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of anti-top quark.')]
            daughter_t1_W = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], passed[i]), PID.w_plus, MODEL) for i in tqdm.trange(len(passed), desc='Finding daughters of W+ boson.')]
            daughter_t2_W = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], passed[i]), PID.w_minus, MODEL) for i in tqdm.trange(len(passed), desc='Finding daughters of W- boson.')]
            daughter_t1_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in daughter_t1_W])
            daughter_t2_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in daughter_t2_W])

            daughter_t1_W_1 = daughter_t1_W_idx[:,0]
            daughter_t1_W_2 = daughter_t1_W_idx[:,1]
            daughter_t1_b = np.array([a["daughter_2_idx"] for a in daughter_t1])

            daughter_t2_W_1 = daughter_t2_W_idx[:,0]
            daughter_t2_W_2 = daughter_t2_W_idx[:,1]
            daughter_t2_b = np.array([a["daughter_1_idx"] for a in daughter_t2])
            
            daughter_h = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], i), PID.higgs, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of higgs boson.')]
            daughter_h_b_1 = np.array([a["daughter_1_idx"] for a in daughter_h])
            daughter_h_b_2 = np.array([a["daughter_2_idx"] for a in daughter_h])
        elif MODEL == 'four_top':
            daughter_t = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], i), PID.top, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of top quark')]
            daughter_anti_t = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], i), PID.anti_top, MODEL) for i in tqdm.tqdm(passed, total=len(passed), desc='Finding daughters of anti-top quark')]
            daughter_t_W = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], passed[i]), PID.w_plus, MODEL) for i in tqdm.trange(len(passed), desc='Finding daughters of W+ boson')]
            daughter_anti_t_W = [process_methods.daughter_finder(helper.to_dataframe(dataset["particle"], passed[i]), PID.w_minus, MODEL) for i in tqdm.trange(len(passed), desc='Finding daughters of W- boson')]
            daughter_t_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_1_idx", "daughter_1_2_idx", "daughter_2_1_idx", "daughter_2_2_idx"]] for dic in daughter_t_W])
            daughter_anti_t_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_1_idx", "daughter_1_2_idx", "daughter_2_1_idx", "daughter_2_2_idx"]] for dic in daughter_anti_t_W])

            daughter_t1_W_1 = daughter_t_W_idx[:,0]
            daughter_t1_W_2 = daughter_t_W_idx[:,1]
            daughter_t2_W_1 = daughter_t_W_idx[:,2]
            daughter_t2_W_2 = daughter_t_W_idx[:,3]
            daughter_t1_b = np.array([a["daughter_1_2_idx"] for a in daughter_t])
            daughter_t2_b = np.array([a["daughter_2_2_idx"] for a in daughter_t])

            daughter_t3_W_1 = daughter_anti_t_W_idx[:,0]
            daughter_t3_W_2 = daughter_anti_t_W_idx[:,1]
            daughter_t4_W_1 = daughter_anti_t_W_idx[:,2]
            daughter_t4_W_2 = daughter_anti_t_W_idx[:,3]
            daughter_t3_b = np.array([a["daughter_1_2_idx"] for a in daughter_anti_t])
            daughter_t4_b = np.array([a["daughter_2_2_idx"] for a in daughter_anti_t])
    else: ### MULTIPROCESSING not yet implemented for ttbar_lep
        if MODEL == 'ttbar' or MODEL == 'ttbar_lep_left' or MODEL =='ttbar_lep_right':
            mp_src_top_1 = [[helper.to_dataframe(dataset["particle"], i), PID.top, MODEL] for i in tqdm.tqdm(passed, total=len(passed), desc='preparing data for multiprocessing.(1/4)')]
            mp_src_top_2 = [[helper.to_dataframe(dataset["particle"], i), PID.anti_top, MODEL] for i in tqdm.tqdm(passed, total=len(passed), desc='preparing data for multiprocessing.(2/4)')]
            mp_src_W_1 = [[helper.to_dataframe(dataset["particle"], passed[i]), PID.w_plus, MODEL] for i in tqdm.trange(len(passed), desc='preparing data for multiprocessing.(3/4)')]
            mp_src_W_2 = [[helper.to_dataframe(dataset["particle"], passed[i]), PID.w_minus, MODEL] for i in tqdm.trange(len(passed), desc='preparing data for multiprocessing.(4/4)')]
            
            process_list = [mp_src_top_1, mp_src_top_2, mp_src_W_1, mp_src_W_2]
            result_name = ["result_t1", "result_t2", "result_t1_W", "result_t2_W"]
            result = OrderedDict()
            print("Start multiprocessing.")
            start = time.time()
            for _process, _name in zip(process_list, result_name):
                print(f"Computing Job: {_name}.")
                with mp.Pool(PROCESS) as p:
                    result[_name] = p.starmap(process_methods.daughter_finder, _process)
                    p.close()
                    p.join()
            print(f"Multiprocessing complete. Cost: {time.time() - start:.3f} s")

            daughter_t1_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in result["result_t1_W"]])
            daughter_t2_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in result["result_t2_W"]])

            daughter_t1_W_1 = daughter_t1_W_idx[:,0]
            daughter_t1_W_2 = daughter_t1_W_idx[:,1]
            daughter_t2_W_1 = daughter_t2_W_idx[:,0]
            daughter_t2_W_2 = daughter_t2_W_idx[:,1]
            daughter_t1_b = np.array([a["daughter_2_idx"] for a in result["result_t1"]])
            daughter_t2_b = np.array([a["daughter_1_idx"] for a in result["result_t2"]])

        elif MODEL == 'ttH':
            mp_src_top_1 = [[helper.to_dataframe(dataset["particle"], i), PID.top, MODEL] for i in tqdm.tqdm(passed, total=len(passed), desc='preparing data for multiprocessing.(1/5)')]
            mp_src_top_2 = [[helper.to_dataframe(dataset["particle"], i), PID.anti_top, MODEL] for i in tqdm.tqdm(passed, total=len(passed), desc='preparing data for multiprocessing.(2/5)')]
            mp_src_W_1 = [[helper.to_dataframe(dataset["particle"], passed[i]), PID.w_plus, MODEL] for i in tqdm.trange(len(passed), desc='preparing data for multiprocessing.(3/5)')]
            mp_src_W_2 = [[helper.to_dataframe(dataset["particle"], passed[i]), PID.w_minus, MODEL] for i in tqdm.trange(len(passed), desc='preparing data for multiprocessing.(4/5)')]
            mp_src_h = [[helper.to_dataframe(dataset["particle"], i), PID.higgs, MODEL] for i in tqdm.tqdm(passed, total=len(passed), desc='preparing data for multiprocessing.(5/5)')]
            
            process_list = [mp_src_top_1, mp_src_top_2, mp_src_W_1, mp_src_W_2, mp_src_h]
            result_name = ["result_t1", "result_t2", "result_t1_W", "result_t2_W", "result_h"]
            result = OrderedDict()
            print("Start multiprocessing.")
            start = time.time()
            for _process, _name in zip(process_list, result_name):
                print(f"Computing Job: {_name}.")
                with mp.Pool(PROCESS) as p:
                    result[_name] = p.starmap(process_methods.daughter_finder, _process)
                    p.close()
                    p.join()
            print(f"Multiprocessing complete. Cost: {time.time() - start:.3f} s")
            daughter_t1_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in result["result_t1_W"]])
            daughter_t2_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_idx", "daughter_2_idx"]] for dic in result["result_t2_W"]])
            daughter_t1_W_1 = daughter_t1_W_idx[:,0]
            daughter_t1_W_2 = daughter_t1_W_idx[:,1]
            daughter_t2_W_1 = daughter_t2_W_idx[:,0]
            daughter_t2_W_2 = daughter_t2_W_idx[:,1]
            daughter_t1_b = np.array([a["daughter_2_idx"] for a in result["result_t1"]])
            daughter_t2_b = np.array([a["daughter_1_idx"] for a in result["result_t2"]])         
            daughter_h_b_1 = np.array([a["daughter_1_idx"] for a in result["result_h"]])
            daughter_h_b_2 = np.array([a["daughter_2_idx"] for a in result["result_h"]])
        elif MODEL == 'four_top':
            mp_src_top = [[helper.to_dataframe(dataset["particle"], i), PID.top, MODEL] for i in tqdm.tqdm(passed, total=len(passed), desc='preparing data for multiprocessing.(1/4)')]
            mp_src_anti_top = [[helper.to_dataframe(dataset["particle"], i), PID.anti_top, MODEL] for i in tqdm.tqdm(passed, total=len(passed), desc='preparing data for multiprocessing.(2/4)')]
            mp_src_top_W = [[helper.to_dataframe(dataset["particle"], passed[i]), PID.w_plus, MODEL] for i in tqdm.trange(len(passed), desc='preparing data for multiprocessing.(3/4)')]
            mp_src_anti_top_W = [[helper.to_dataframe(dataset["particle"], passed[i]), PID.w_minus, MODEL] for i in tqdm.trange(len(passed), desc='preparing data for multiprocessing.(4/4)')]
            
            process_list = [mp_src_top, mp_src_anti_top, mp_src_top_W, mp_src_anti_top_W]
            result_name = ["result_top", "result_anti_top", "result_top_W", "result_anti_top_W"]
            result = OrderedDict()
            print("Start multiprocessing.")
            start = time.time()
            for _process, _name in zip(process_list, result_name):
                print(f"Computing Job: {_name}.")
                with mp.Pool(PROCESS) as p:
                    result[_name] = p.starmap(process_methods.daughter_finder, _process)
                    p.close()
                    p.join()
            print(f"Multiprocessing complete. Cost: {time.time() - start:.3f} s")
            
            daughter_t_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_1_idx", "daughter_1_2_idx", "daughter_2_1_idx", "daughter_2_2_idx"]] for dic in result["result_top_W"]])
            daughter_anti_t_W_idx = np.array([[ dic[item] for item in dic if item in ["daughter_1_1_idx", "daughter_1_2_idx", "daughter_2_1_idx", "daughter_2_2_idx"]] for dic in result["result_anti_top_W"]])
            daughter_t1_W_1 = daughter_t_W_idx[:,0]
            daughter_t1_W_2 = daughter_t_W_idx[:,1]
            daughter_t2_W_1 = daughter_t_W_idx[:,2]
            daughter_t2_W_2 = daughter_t_W_idx[:,3]
            daughter_t3_W_1 = daughter_anti_t_W_idx[:,0]
            daughter_t3_W_2 = daughter_anti_t_W_idx[:,1]
            daughter_t4_W_1 = daughter_anti_t_W_idx[:,2]
            daughter_t4_W_2 = daughter_anti_t_W_idx[:,3]
            daughter_t1_b = np.array([a["daughter_1_2_idx"] for a in result["result_top"]])
            daughter_t2_b = np.array([a["daughter_2_2_idx"] for a in result["result_top"]])
            daughter_t3_b = np.array([a["daughter_1_1_idx"] for a in result["result_anti_top"]])
            daughter_t4_b = np.array([a["daughter_2_1_idx"] for a in result["result_anti_top"]])

    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording the kinematics variables of partons and jets in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    if MODEL == 'ttbar' or MODEL in require_lepton:
        parton_idx = np.stack((daughter_t1_b, daughter_t1_W_1, daughter_t1_W_2, daughter_t2_b, daughter_t2_W_1, daughter_t2_W_2), axis=1)
    elif MODEL == 'ttH':
        parton_idx = np.stack((daughter_t1_b, daughter_t1_W_1, daughter_t1_W_2, daughter_t2_b, daughter_t2_W_1, daughter_t2_W_2, daughter_h_b_1, daughter_h_b_2), axis=1)
    elif MODEL == 'four_top':
        parton_idx = np.stack((daughter_t1_b, daughter_t1_W_1, daughter_t1_W_2, 
                               daughter_t2_b, daughter_t2_W_1, daughter_t2_W_2, 
                               daughter_t3_b, daughter_t3_W_1, daughter_t3_W_2,
                               daughter_t4_b, daughter_t4_W_1, daughter_t4_W_2), axis=1)

    sourece_features = ["PT", "Eta", "Phi", "Mass", "PID"]
    storage_name = ["pt", "eta", "phi", "mass", "pdgid"]

    parton_features = OrderedDict()
    for a, b in tqdm.tqdm(zip(sourece_features, storage_name), total=(len(sourece_features)), desc="Storing parton's kinematics information"):
        parton_features[b] = np.array([helper.fetch_kinematics_properties_from_dataset(helper.to_dataframe(dataset["particle"], passed[i]), parton_idx[i], a) for i in range(len(passed))])

    parton_barcode = np.tile(barcode, (len(passed),1))
    parton_features["parton_barcode"] = parton_barcode

    for a in dataset.keys():
        for b in dataset[a].keys():
            dataset[a][b] = dataset[a][b][passed]
    if MODEL == 'ttbar_lep' or MODEL == 'ttbar_lep_left':
        parton_masks = parton_barcode != 40
    elif MODEL == "ttbar_lep_right":
        parton_masks = parton_barcode != 20
    else :
        parton_masks = parton_barcode != 0

    parton_features['masks'] = parton_masks
    
    jet_features_list = ['event', 'pt', 'eta', 'phi', 'btag', 'mass', 'num_of_jets', 'charge']
    jet_features = OrderedDict()
    for feature in tqdm.tqdm(jet_features_list, total=len(jet_features_list), desc='Storing jet informations.'):
        if feature not in ['event', 'num_of_jets'] :
            jet_features[feature] = np.array([np.pad(np.array(x, dtype=np.float64), (0, MAX_NUM_OF_JETS - len(x)), 'constant', constant_values=(0, -999)).astype(np.float64).tolist() for x in dataset['jet'][feature]])
        else:
            jet_features[feature] = dataset['jet'][feature]
    print("+------------------------------------------------------------------------------------------------------+")
    print("Finished to record the kinematics variables of partons and jets in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    
    print("+------------------------------------------------------------------------------------------------------+")
    print("Starting parton-jet matching.")
    print("+------------------------------------------------------------------------------------------------------+")
    truth_matching_result = [deltaR_matching(NUM_OF_PARTON, 
                      dataset['jet']['num_of_jets'][i], 
                      parton_features["eta"][i], 
                      parton_features["phi"][i], 
                      dataset['jet']['eta'][i], 
                      dataset['jet']['phi'][i],
                      0.4, 
                      MODEL) for i in tqdm.trange(len(passed), desc='Computing truth matching')]
    jet_parton_index = np.array([np.pad(x[0], (0, MAX_NUM_OF_JETS - len(x[0])), 'constant', constant_values=(0, -999)).tolist() for x in truth_matching_result])
    parton_jet_index = np.array([x[1].tolist() for x in truth_matching_result])
    print("+------------------------------------------------------------------------------------------------------+")
    print("Parton-jet matching finished.")
#     print("Parton-jet matching finished. Cost: {0:.1f} s".format(time.time()-start))
    print("+------------------------------------------------------------------------------------------------------+")
       
    if MODEL == 'four_top' and COMPUTE_CHI2:
        while True:
            print("Computing chi-square for four top model might be very slow, please check before continue.")
            ans = input("Do you want to comput chi-square for four top model.(y,n)")
            if ans in ['y', 'Y']:
                break
            elif ans in ['n', 'N']:
                COMPUTE_CHI2 = False
                break
            else:
                print("Invalid answer.")
    if MODEL in require_lepton and COMPUTE_CHI2:
        print("Chi-square minimizer is not yet support semi-leptonic channel. This step will be skipped directly.")
        COMPUTE_CHI2 = False
    if COMPUTE_CHI2:
        print("+------------------------------------------------------------------------------------------------------+")
        print("Starting chi-square matching.")
        print("+------------------------------------------------------------------------------------------------------+")
                
        if int(PROCESS) == 1:
            EXTRA = kargs["EXTRA"]
            chi_value = []
            parton_jet_index_chi2 = []
            jet_parton_index_chi2 = []
            least_10_chi2_candidate = []
            least_10_chi2_value = []

            for i in tqdm.trange(len(passed), desc="Computing chi-square minimizer"):
                tmp_chi2_result, tmp_parton_jet_index, tmp_jet_parton_index, tmp_least_10_chi2_candidate, tmp_least_10_chi2_value = chi_square_minimizer(jet_features['pt'][i], 
                                                                                                                                     jet_features['eta'][i], 
                                                                                                                                     jet_features['phi'][i], 
                                                                                                                                     jet_features['btag'][i], 
                                                                                                                                     jet_features['mass'][i], 
                                                                                                                                     MODEL, 
                                                                                                                                     EXTRA, 
                                                                                                                                     jet_features['num_of_jets'][i])
                chi_value.append(tmp_chi2_result)
                parton_jet_index_chi2.append(tmp_parton_jet_index)
                jet_parton_index_chi2.append(tmp_jet_parton_index)    
                least_10_chi2_candidate.append(tmp_least_10_chi2_candidate)    
                least_10_chi2_value.append(tmp_least_10_chi2_value)    
            chi_value = np.array(chi_value)
            parton_jet_index_chi2 = np.array(parton_jet_index_chi2)
            jet_parton_index_chi2 = np.array(jet_parton_index_chi2)
            least_10_chi2_candidate = np.array(least_10_chi2_candidate)
            least_10_chi2_value = np.array(least_10_chi2_value)
            if (chi_value == -1).sum() != 0: print("There exist some events failed to compute chi-square reconstuction.")
        else: 
            EXTRA = kargs["EXTRA"]
            src_chi2 = [[jet_features['pt'][i], 
                         jet_features['eta'][i], 
                         jet_features['phi'][i], 
                         jet_features['btag'][i], 
                         jet_features['mass'][i], 
                         MODEL, 
                         EXTRA, 
                         jet_features['num_of_jets'][i]] for i in tqdm.trange(len(passed), desc="Preparing data for multiprocessing.")]
            print("Start multiprocessing.")
            start = time.time()
            with mp.Pool(PROCESS) as p:
                _result_chi2 = p.starmap(chi_square_minimizer, src_chi2)
                p.close()
                p.join()
            print(f"Multiprocessing complete. Cost: {time.time() - start:.3f} s")
            _result_chi2 = np.array(_result_chi2, dtype='object')
            
            chi_value = np.array(_result_chi2[:, 0], dtype=np.float64)
            parton_jet_index_chi2 = np.array([ x for x in _result_chi2[:, 1]])
            jet_parton_index_chi2 = np.array([ np.pad(x, (0, MAX_NUM_OF_JETS - len(x)), 'constant', constant_values=(0, -999)).tolist() for x in _result_chi2[:, 2]])
            least_10_chi2_candidate = np.array([ x for x in _result_chi2[:, 3]])
            least_10_chi2_value = np.array([ x for x in _result_chi2[:, 4]])
            if (chi_value == -1).sum() != 0: print("There exist some events failed to compute chi-square reconstuction.")
        print("+------------------------------------------------------------------------------------------------------+")
        print("Chi-square matching finished.")
        print("+------------------------------------------------------------------------------------------------------+")
    if MODEL == 'ttbar_lep' or MODEL == 'ttbar_lep_left' or MODEL == "ttbar_lep_right":
        lepton_pt = np.array([float(a) if len(a) != 0 else float(b) for a, b in zip(dataset['muon']['pt'], dataset['electron']['pt'])])
        lepton_eta = np.array([float(a) if len(a) != 0 else float(b) for a, b in zip(dataset['muon']['eta'], dataset['electron']['eta'])])
        lepton_phi = np.array([float(a) if len(a) != 0 else float(b) for a, b in zip(dataset['muon']['phi'], dataset['electron']['phi'])])
        lepton_features = OrderedDict((
            ("pt", lepton_pt),
            ("eta", lepton_eta),
            ("phi", lepton_phi),
        ))
        MET_features = OrderedDict((
            ("MET", np.array([x for x in dataset['MissingET']['MET']])),
            ("eta", np.array([x for x in dataset['MissingET']['eta']])),
            ("phi", np.array([x for x in dataset['MissingET']['phi']])),
        ))
        print("+------------------------------------------------------------------------------------------------------+")
        print("Starting lepton and neutrino matching.")
        print("+------------------------------------------------------------------------------------------------------+")
        if MODEL == 'ttbar_lep' or MODEL == 'ttbar_lep_left':
            lepton_delta_R_result = np.zeros(len(passed))
            met_delta_phi_result = np.zeros(len(passed))
            for i in tqdm.trange(len(passed), desc='Computing delta R for leptons'):
                _delta_R_lepton = delta_R(parton_features['eta'][i][1], parton_features['phi'][i][1], lepton_features['eta'][i], lepton_features['phi'][i])
                _delta_phi_met = deltaPhi(parton_features['phi'][i][2], MET_features['phi'][i])
                if _delta_R_lepton < 0.4:
                    lepton_delta_R_result[i] = 1
                else : 
                    lepton_delta_R_result[i] = 0
                if _delta_phi_met < 0.4:
                    met_delta_phi_result[i] = 1
                else : 
                    met_delta_phi_result[i] = 0
        elif MODEL == 'ttbar_lep_right':
            lepton_delta_R_result = np.zeros(len(passed))
            met_delta_phi_result = np.zeros(len(passed))
            for i in tqdm.trange(len(passed), desc='Computing delta R for leptons'):
                _delta_R_lepton = delta_R(parton_features['eta'][i][5], parton_features['phi'][i][4], lepton_features['eta'][i], lepton_features['phi'][i])
                _delta_phi_met = deltaPhi(parton_features['phi'][i][4], MET_features['phi'][i])
                if _delta_R_lepton < 0.4:
                    lepton_delta_R_result[i] = 1
                else : 
                    lepton_delta_R_result[i] = 0
                if _delta_phi_met < 0.4:
                    met_delta_phi_result[i] = 1
                else : 
                    met_delta_phi_result[i] = 0
        print("+------------------------------------------------------------------------------------------------------+")
        print("Neutrino and lepton matching finished.")
        print("+------------------------------------------------------------------------------------------------------+") 
    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording barcode information.")
    print("+------------------------------------------------------------------------------------------------------+")
    jet_barcode = jet_parton_index.copy()
    for i in range(len(barcode)):
        jet_barcode = np.where(jet_barcode == i, barcode[i], jet_barcode) 
    jet_features['barcode'] = jet_barcode
    if COMPUTE_CHI2:
        jet_barcode_chi2 = jet_parton_index_chi2.copy()
        for i in range(len(barcode)):
            jet_barcode_chi2 = np.where(jet_barcode_chi2 == i, barcode[i], jet_barcode_chi2) 
        jet_features['barcode'] = jet_barcode_chi2
    print("+------------------------------------------------------------------------------------------------------+")
    print("Barcode information has beed record.")
    print("+------------------------------------------------------------------------------------------------------+")

    print("+------------------------------------------------------------------------------------------------------+")
    print("Writing event record to the hdf5 file.")
    print("+------------------------------------------------------------------------------------------------------+")
    if OUTPUT_FILE.split(".")[-1] == "h5":
        pass
    else:
        OUTPUT_FILE = OUTPUT_FILE+".h5"
    if MODEL == 'ttbar' or MODEL in require_lepton:
        targets = OrderedDict((
            ("left_target", parton_jet_index.T[:3]),
            ("right_target", parton_jet_index.T[3:]),
        ))
        masks = OrderedDict((
            ("left_target",  np.any(parton_jet_index[:,:3] < 0, 1)),
            ("right_target",  np.any(parton_jet_index[:,3:] < 0, 1)),
        ))
    elif MODEL == 'ttH':
        targets = OrderedDict((
            ("left_target", parton_jet_index.T[:3]),
            ("right_target", parton_jet_index.T[3:6]),
            ("higgs_target", parton_jet_index.T[6:]),
        ))
        masks = OrderedDict((
            ("left_target",  np.any(parton_jet_index[:,:3] < 0, 1)),
            ("right_target",  np.any(parton_jet_index[:,3:6] < 0, 1)),
            ("higgs_target",  np.any(parton_jet_index[:,6:] < 0, 1)),
        ))
    elif MODEL == 'four_top':
        targets = OrderedDict((
            ("first_target", parton_jet_index.T[:3]),
            ("second_target", parton_jet_index.T[3:6]),
            ("third_target", parton_jet_index.T[6:9]),
            ("fourth_target", parton_jet_index.T[9:12]),
        ))
        masks = OrderedDict((
            ("first_target",  np.any(parton_jet_index[:,:3] < 0, 1)),
            ("second_target",  np.any(parton_jet_index[:,3:6] < 0, 1)),
            ("third_target",  np.any(parton_jet_index[:,6:9] < 0, 1)),
            ("fourth_target",  np.any(parton_jet_index[:,9:12] < 0, 1)),
        ))
    if COMPUTE_CHI2:
        chi_square_result = OrderedDict((
            ("chi_value", chi_value),
            ("least_10_chi2_candidate", least_10_chi2_candidate),
            ("least_10_chi2_value", least_10_chi2_value),
        ))
    if MODEL in require_lepton: 
        output = IO_module(OUTPUT_FILE, 
            MODEL, 
            targets=targets, 
            masks=masks, 
            parton_features=parton_features, 
            jet_features=jet_features, 
            MET_features=MET_features,
            lepton_features=lepton_features,
            usage='parse',
            contain_chi_square=False,
           )
    else:
        if COMPUTE_CHI2:
            output = IO_module(OUTPUT_FILE, 
                MODEL, 
                targets=targets, 
                masks=masks, 
                parton_features=parton_features, 
                jet_features=jet_features, 
                chi_square_result=chi_square_result,
                usage='parse',
                contain_chi_square=True,
               )
        else:
            output = IO_module(OUTPUT_FILE, 
                MODEL, 
                targets=targets, 
                masks=masks, 
                parton_features=parton_features, 
                jet_features=jet_features, 
                usage='parse',
                contain_chi_square=False,
               )
    output.wrtite_hdf5()
    if os.path.isfile(OUTPUT_FILE):
        print("+------------------------------------------------------------------------------------------------------+")
        print("Event record has been send to {0}.".format(OUTPUT_FILE))
        print("+------------------------------------------------------------------------------------------------------+")
    else:
        print("+------------------------------------------------------------------------------------------------------+")
        print("Error occur. Fail to write event record to {0}.".format(OUTPUT_FILE))
        print("+------------------------------------------------------------------------------------------------------+")
