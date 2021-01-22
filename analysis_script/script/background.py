"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
#Import packages
import uproot, time
import pandas as pd 
import numpy as np 
from .particle_properties_uproot import particle_properties  #import particle properties helper function from particle_properties.py
from .jet_properties_uproot import jet_properties  #import jet properties helper function from jet_properties.py
import h5py, sys, traceback, os, tqdm
from .utilize import delta_R, deltaPhi, pdgid, event_selection, quark_finder, particle_tracing, deltaR_matching, barcode_recorder, chi_square_minimizer
import multiprocessing as mp

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


def background(INPUT, OUTPUT, MODEL, SINGLE, PROCESS):
    PID = pdgid()

    if int(SINGLE) == 1:
        try:
            print("Loading root file from {0}.".format(INPUT_FILE))
            data  = uproot.open(INPUT_FILE)['Delphes']
        except:
            print('Please check input file path.')
        particle = particle_properties(data)
        jet = jet_properties(data)

    elif int(SINGLE) != 1 and bool(SINGLE.isdigit) == True:
        _data = []
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

        particle = properties_for_particle(_particle_event, _particle_pt, _particle_eta, _particle_phi, _particle_pid, _particle_M1, _particle_M2, _particle_D1, _particle_D2, _particle_status, _particle_rapidity, _particle_mass, _particle_charge)

        jet = properties_for_jet(_jet_event, _jet_pt, _jet_eta, _jet_phi, _jet_btag, _jet_area, _jet_mass, _jet_charge)
    else:
        print("Please inpurt a correct mode.\n1.-s 1\n2.-s > 1")
    
    barcode = np.array([34, 40, 40, 17, 20, 20])
    print("+------------------------------------------------------------------------------------------------------+")
    print("Start jet selection.")
    print("+------------------------------------------------------------------------------------------------------+")
    marker_event, marker_jet, marker_btag = event_selection(jet.pt, jet.eta, jet.btag, MODEL)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Jet selection doen. {0} events has been selected.".format(np.sum(marker_event == 1)))
    print("+------------------------------------------------------------------------------------------------------+")
    
    print("+------------------------------------------------------------------------------------------------------+")
    print("Recording the kinematics variables of jets in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    #Record the kinematical variables of jet in the selected event.
    jet_pt = []
    jet_eta = []
    jet_phi = []
    jet_mass = []
    jet_btag = []

    for i in tqdm.trange(len(jet.event)):
        if marker_event[i] == 1:
            jet_pt_tmp = []
            jet_eta_tmp = []
            jet_phi_tmp = []
            jet_mass_tmp = []
            jet_btag_tmp = []
            for j in range(len(jet.pt[i])):
                jet_pt_tmp.append(jet.pt[i][j])
                jet_eta_tmp.append(jet.eta[i][j])
                jet_phi_tmp.append(jet.phi[i][j])
                jet_mass_tmp.append(jet.mass[i][j])
                jet_btag_tmp.append(jet.btag[i][j])
                    
            jet_pt.append(jet_pt_tmp)
            jet_eta.append(jet_eta_tmp)
            jet_phi.append(jet_phi_tmp)
            jet_mass.append(jet_mass_tmp)
            jet_btag.append(jet_btag_tmp)
    print("+------------------------------------------------------------------------------------------------------+")
    print("Finished to record the kinematics variables of jets in the selected event.")
    print("+------------------------------------------------------------------------------------------------------+")
    
    print("+------------------------------------------------------------------------------------------------------+")
    print("Starting chi-square matching.")
    print("+------------------------------------------------------------------------------------------------------+")
    start = time.time()
    chi2_value = []
    chi2_result = []
    _src_chi2 = []
    for i in range(len(jet_pt)):
        _src_chi2.append([jet_pt[i], jet_eta[i], jet_phi[i], jet_btag[i], jet_mass[i], MODEL, EXTRA])
    print("Using {0} process for accelerating speed.".format(PROCESS))
    with mp.Pool(PROCESS) as p:
        _result_chi2 = p.starmap(chi_square_minimizer, _src_chi2)
        p.close()
        p.join()
    for i in range(len(_result_chi2)):
        if _result_chi2[i][0] == -1 :
            print("+++++++++++++++++++++++++++++++++++++")
            print('An error occur!')
            print("Ebvent number: {0}, chi2 vslue: {1}, _tmp_parton_jet_index: {2}, _ttmp_jet_parton_index: {3}".format(i, _result_chi2[i][0], _result_chi2[i][1], _result_chi2[i][2]))
            print("+++++++++++++++++++++++++++++++++++++")
        if _result_chi2[i][0] != -1 :
            chi2_value.append(_result_chi2[i][0])
            chi2_result.append(_result_chi2[i][1])

    print("+------------------------------------------------------------------------------------------------------+")
    print("Chi-square matching finished. Cost: {0:.1f} s".format(time.time() - start))
    print("+------------------------------------------------------------------------------------------------------+")
    
    # Save selected events
    lene = len(parton_pdgid)

    print("+------------------------------------------------------------------------------------------------------+")
    print("Writing event record to the hdf5 file.")
    print("+------------------------------------------------------------------------------------------------------+")

    # Save the event which pass the selection
    with h5py.File(OUTPUT_FILE,'w') as f:
        dt = h5py.vlen_dtype(np.dtype('float16'))
        
        hdf5_chi2_result = f.create_dataset('chi2_result', (lene, ), dtype=dt)
        hdf5_jet_pt = f.create_dataset('jet_pt', (lene, ), dtype=dt)
        hdf5_jet_eta = f.create_dataset('jet_eta', (lene, ), dtype=dt)
        hdf5_jet_phi = f.create_dataset('jet_phi', (lene, ), dtype=dt)
        hdf5_jet_mass = f.create_dataset('jet_mass', (lene, ), dtype=dt)
        hdf5_jet_btag = f.create_dataset('jet_btag', (lene, ), dtype=dt)
        print("Writing jet information to {0}".format(OUTPUT_FILE))
        
        for i in tqdm.trange(lene):
            hdf5_chi2_result[i] = chi2_result[i]
            hdf5_jet_pt[i] = jet_pt[i]
            hdf5_jet_eta[i] = jet_eta[i]
            hdf5_jet_phi[i] = jet_phi[i]
            hdf5_jet_mass[i] = jet_mass[i]
            hdf5_jet_btag[i] = jet_btag[i]

    print("+------------------------------------------------------------------------------------------------------+")
    print("Event record has been send to {0}.".format(OUTPUT_FILE))
    print("+------------------------------------------------------------------------------------------------------+")