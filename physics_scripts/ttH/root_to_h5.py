#!/usr/bin/env python3

from ssl import HAS_SNI
from symbol import star_expr
import awkward as ak
import h5py
import numpy as np
import tqdm
import uproot


DATA_STRUCTURE = np.dtype(
    [
        ("bhad", "int32"),
        ("blep", "int32"),
        ("q1", "int32"),
        ("q2", "int32"),
        ("b1", "int32"),
        ("b2", "int32"),
        ("lepton_pt1", "f4"),
        ("lepton_phi1", "f4"),
        ("lepton_eta1", "f4"),
        ("lepton_cl_eta1", "f4"),
        ("lepton_m1", "f4"),
        ("lepton_charge1", "f4"),
        ("lepton_pid1", "int32"),
        ("met_met1", "f4"),
        ("met_phi1", "f4"),
        ("sumet1", "f4"),
        ("Bhad", "int32"),
        ("Blep", "int32"),
        ("Q1", "int32"),
        ("Q2", "int32"),
        ("B1", "int32"),
        ("B2", "int32"),
    ]
)

DATA_STRUCTUREA = np.dtype(
    [
        ("higgsm", "f4"),
        ("Ravgbb", "f4"),
        ("Rmaxptbb", "f4"),
        ("etamaxjj", "f4"),
        ("massminRbb", "f4"),
        ("massminRjj", "f4"),
        ("Nbbhiggs30", "f4"),
        ("HThad", "f4"),
        ("Rminlbb", "f4"),
        ("mhblep", "f4"),
        ("Rhbb", "f4"),
        ("Rhtt", "f4"),
        ("Rhtlep", "f4"),
        ("Rhbhad", "f4"),
        ("likelihoodmax", "f4"),
    ]
)

JET_STRUCTURE = np.dtype(
    [
        ("jet_pt1", "float64"),
        ("jet_phi1", "float64"),
        ("jet_eta1", "float64"),
        ("jet_m1", "float64"),
        ("jet_barcode1", "float64"),
        ("jet_btag1", "float64"),
    ]
)

PARTON_STRUCTURE = np.dtype(
    [
        ("parton_pt1", "float64"),
        ("parton_phi1", "float64"),
        ("parton_eta1", "float64"),
        ("parton_m1", "float64"),
        ("parton_pid1", "int32"),
        ("parton_barcode1", "int32"),
    ]
)
def get_flattened_eventA(event):
    """Returns event in the specified data structure."""
    flattened_event = np.zeros(1, dtype=DATA_STRUCTUREA)

   

    flattened_event["higgsm"] = event["higgsm"]
    flattened_event["Ravgbb"] = event["Ravgbb"]
    flattened_event["Rmaxptbb"] = event["Rmaxptbb"]
    flattened_event["etamaxjj"] = event["etamaxjj"]
    flattened_event["massminRbb"] = event["massminRbb"]
    flattened_event["massminRjj"] = event["massminRjj"]
    flattened_event["Nbbhiggs30"] = event["Nbbhiggs30"]
    flattened_event["HThad"] = event["HThad"]
    flattened_event["Rminlbb"] = event["Rminlbb"]
    flattened_event["mhblep"] = event["mhblep"]
    flattened_event["Rhbb"] = event["Rhbb"]
    flattened_event["Rhtt"] = event["Rhtt"]
    flattened_event["Rhtlep"] = event["Rhtlep"]
    flattened_event["Rhbhad"] = event["Rhbhad"]
    flattened_event["likelihoodmax"] = event["likelihoodmax"]

    return flattened_event


def create_structured_dataA(file_list: str, tree_name: str, n_evt_max: int = None):
    """Reads tree from a list of files and creates structured data.

    This function takes a list of files, concatenates the trees of the specified
    name, applies an event selection according to the passes_selection() method
    and then creates a flattened data structure according to the
    get_flattened_event() method. Returns a numpy array of all selected and
    flattened events. n_evt_max can be used to restrict the number of processed
    events.

    """
    tree = uproot.concatenate(
        f"{file_list}:{tree_name}",
        [
            "higgsm",
            "Ravgbb",
            "Rmaxptbb",
            "etamaxjj",
            "massminRbb",
            "massminRjj",
            "Nbbhiggs30",
            "HThad",
            "Rminlbb",
            "mhblep",
            "Rhbb",
            "Rhtt",
            "Rhtlep",
            "Rhbhad",
            "likelihoodmax",
        ],
    )

    n_evt_tot = min(int(n_evt_max), len(tree)) if n_evt_max else len(tree)
    n_evt_sel = 0

    data = np.zeros(n_evt_tot, dtype=DATA_STRUCTUREA)

    print(
        f'Reading {n_evt_tot:,} events from concatenated '
        f'tree "{tree_name}" with {len(tree):,} events.'
    )

    for n_evt in tqdm.tqdm(range(n_evt_tot)):
        event = tree[n_evt]
        #if not passes_selection(event):
        #    continue
        data[n_evt_sel] = get_flattened_eventA(event)
        n_evt_sel += 1

    print(
        f"Selected {n_evt_sel:,} out of {n_evt_tot:,} events "
        f"({n_evt_sel / n_evt_tot * 100:.1f}%%)."
    )

    return data[:n_evt_sel]

def passes_selection(event):
    """Checks whether an event passes the selection criteria.

    The function selects events with at least 3 jets, exactly 3 leptons and at
    least one jet tagged at the 77% DL1r operating point.

    """
    if not len(event["jet_pt"]) >= 3:
        return False

    if not len(event["Leptons_Pt"]) == 3:
        return False

    if not ak.sum(event["jet_tagWeightBin_DL1r_Continuous"] >= 2) >= 1:
        return False

    return True


def get_flattened_event(event):
    """Returns event in the specified data structure."""
    flattened_event = np.zeros(1, dtype=DATA_STRUCTURE)
    jet_event = np.zeros((1,18), dtype=JET_STRUCTURE)
    parton_event = np.zeros((1,8), dtype=PARTON_STRUCTURE)
   

    flattened_event["bhad"] = event["bhad"]
    flattened_event["blep"] = event["blep"]
    flattened_event["b1"] = event["b1"]
    flattened_event["b2"] = event["b2"]
    flattened_event["q1"] = event["q1"]
    flattened_event["q2"] = event["q2"]

    flattened_event["Bhad"] = event["target_Bhad"]
    flattened_event["Blep"] = event["target_Blep"]
    flattened_event["B1"] = event["target_Higgs_b1"]
    flattened_event["B2"] = event["target_Higgs_b2"]
    flattened_event["Q1"] = event["target_Q1"]
    flattened_event["Q2"] = event["target_Q2"]

    flattened_event["lepton_pt1"] = event["lepton_pt1"]
    flattened_event["lepton_eta1"] = event["lepton_eta1"]
    flattened_event["lepton_cl_eta1"] = event["lepton_cl_eta1"]
    flattened_event["lepton_phi1"] = event["lepton_phi1"]
    flattened_event["lepton_m1"] = event["lepton_m1"]
    flattened_event["lepton_charge1"] = event["lepton_charge1"]
    flattened_event["lepton_pid1"] = event["lepton_pid1"]

    flattened_event["met_met1"] = event["met_met1"]
    flattened_event["met_phi1"] = event["met_phi1"]
    flattened_event["sumet1"] = event["sumet1"]
    jet_event["jet_phi1"] = event["jet_phi1"]
    jet_event["jet_pt1"] = event["jet_pt1"]
    jet_event["jet_eta1"] = event["jet_eta1"]
    jet_event["jet_m1"] = event["jet_m1"]
    jet_event["jet_barcode1"] = event["jet_barcode1"]
    jet_event["jet_btag1"] = event["jet_has_btag1"]
    parton_event["parton_phi1"] = event["parton_phi1"]
    parton_event["parton_pt1"] = event["parton_pt1"]
    parton_event["parton_pid1"] = event["parton_pid1"]
    parton_event["parton_barcode1"] = event["parton_barcode1"]
    parton_event["parton_m1"] = event["parton_m1"]
    parton_event["parton_eta1"] = event["parton_eta1"]
    return flattened_event, jet_event, parton_event


def create_structured_data(file_list: str, tree_name: str, n_evt_max: int = None):
    """Reads tree from a list of files and creates structured data.

    This function takes a list of files, concatenates the trees of the specified
    name, applies an event selection according to the passes_selection() method
    and then creates a flattened data structure according to the
    get_flattened_event() method. Returns a numpy array of all selected and
    flattened events. n_evt_max can be used to restrict the number of processed
    events.

    """
    tree = uproot.concatenate(
        f"{file_list}:{tree_name}",
        [
            "bhad",
            "blep",
            "q1",
            "q2",
            "b1",
            "b2",
            "lepton_pt1",
            "lepton_eta1",
            "lepton_cl_eta1",
            "lepton_phi1",
            "lepton_m1",
            "lepton_charge1",
            "lepton_pid1",
            "met_met1",
            "met_phi1",
            "sumet1",
            "jet_phi1",
            "jet_pt1",
            "jet_has_btag1",
            "jet_eta1",
            "jet_barcode1",
            "jet_m1",
            "parton_pt1",
            "parton_phi1",
            "parton_pid1",
            "parton_eta1",
            "parton_barcode1",
            "parton_m1",
            "target_Higgs_b1",
            "target_Higgs_b2",
            "target_Blep",
            "target_Bhad",
            "target_Q1",
            "target_Q2",

        ],
    )

    n_evt_tot = min(int(n_evt_max), len(tree)) if n_evt_max else len(tree)
    n_evt_sel = 0

    data = np.zeros(n_evt_tot, dtype=DATA_STRUCTURE)
    jet_data = np.zeros((n_evt_tot, 18), dtype=JET_STRUCTURE)
    parton_data = np.zeros((n_evt_tot, 8), dtype=PARTON_STRUCTURE)

    print(
        f'Reading {n_evt_tot:,} events from concatenated '
        f'tree "{tree_name}" with {len(tree):,} events.'
    )

    for n_evt in tqdm.tqdm(range(n_evt_tot)):
        event = tree[n_evt]
        #if not passes_selection(event):
        #    continue
        data[n_evt_sel], jet_data[n_evt_sel], parton_data[n_evt_sel] = get_flattened_event(event)
        n_evt_sel += 1

    print(
        f"Selected {n_evt_sel:,} out of {n_evt_tot:,} events "
        f"({n_evt_sel / n_evt_tot * 100:.1f}%%)."
    )

    return data[:n_evt_sel], jet_data[:n_evt_sel], parton_data[:n_evt_sel]


def main():
    """Reads signal and background files and converts them to HDF5."""
   
    signal_filesA = f"BDTana.root"
    signal_dataA = create_structured_dataA(signal_filesA, "nominal", 1e7)
    signal_files = f"BDTFormat.root"
    signal_data, jet_data, parton_data = create_structured_data(signal_files, "KLFitter_output", 1e7)

    print(f'Created structured data from input files "{signal_files}".')

    with h5py.File('tth_with_spanet_no_signal.h5', 'r') as fileA:
        h_assignment_probability = fileA["spanet"]["higgs_target"]["assignment_probability"][:]
        h_detection_probability = fileA["spanet"]["higgs_target"]["detection_probability"][:]
        h_marginal_probability = fileA["spanet"]["higgs_target"]["marginal_probability"][:]
        h_b1 = fileA["spanet"]["higgs_target"]["b1"][:]
        h_b2 = fileA["spanet"]["higgs_target"]["b2"][:]
        h_mask = fileA["spanet"]["higgs_target"]["mask"][:]

        l_assignment_probability = fileA["spanet"]["left_target"]["assignment_probability"][:]
        l_detection_probability = fileA["spanet"]["left_target"]["detection_probability"][:]
        l_marginal_probability = fileA["spanet"]["left_target"]["marginal_probability"][:]
        l_b = fileA["spanet"]["left_target"]["b"][:]
        l_mask = fileA["spanet"]["left_target"]["mask"][:]

        r_assignment_probability = fileA["spanet"]["right_target"]["assignment_probability"][:]
        r_detection_probability = fileA["spanet"]["right_target"]["detection_probability"][:]
        r_marginal_probability = fileA["spanet"]["right_target"]["marginal_probability"][:]
        r_b = fileA["spanet"]["right_target"]["b"][:]
        r_q1 = fileA["spanet"]["right_target"]["q1"][:]
        r_q2 = fileA["spanet"]["right_target"]["q2"][:]
        r_mask = fileA["spanet"]["right_target"]["mask"][:]

        s_prediction = fileA["spanet"]["signal"]["prediction"][:]
        s_probability = fileA["spanet"]["signal"]["probability"][:]
        s_target = fileA["spanet"]["signal"]["target"][:]

   
    
    with h5py.File("BDT_KLFitter_output.h5", "w") as file:
        g1=file.create_group("klfitter")
        c1=g1.create_group("left_target")
        c2=g1.create_group("right_target")
        c3=g1.create_group("higgs_target")
        g7=file.create_group("spanet")
        g71 = g7.create_group("signal")
        g71["prediction"] = s_prediction
        g71["probability"] = s_probability
        g71["target"] = s_target

        g72 = g7.create_group("left_target")
        g72["assignment_probability"] = l_assignment_probability
        g72["detection_probability"] =l_detection_probability
        g72["marginal_probability"] = l_marginal_probability
        g72["b"] = l_b
        g72["mask"] = l_mask

        g73 = g7.create_group("right_target")
        g73["assignment_probability"] = r_assignment_probability
        g73["detection_probability"] =r_detection_probability
        g73["marginal_probability"] = r_marginal_probability
        g73["b"] = r_b
        g73["q1"] = r_q1
        g73["q2"] = r_q2
        g73["mask"] = r_mask

        g74 = g7.create_group("higgs_target")
        g74["assignment_probability"] = h_assignment_probability
        g74["detection_probability"] =h_detection_probability
        g74["marginal_probability"] = h_marginal_probability
        g74["b"] = h_b1
        g74["b2"] = h_b2
        g74["mask"] = h_mask

        c2["b"]=signal_data["bhad"]
        c1["b"]=signal_data["blep"]
        c2["q1"]=signal_data["q1"]
        c2["q2"]=signal_data["q2"]
        c3["b1"]=signal_data["b1"]
        c3["b2"]=signal_data["b2"]

        dg1=file.create_group("target")
        dc1=dg1.create_group("left_target")
        dc2=dg1.create_group("right_target")
        dc3=dg1.create_group("higgs_target")
        dc2["b"]=signal_data["Bhad"]
        dc1["b"]=signal_data["Blep"]
        dc2["q1"]=signal_data["Q1"]
        dc2["q2"]=signal_data["Q2"]
        dc3["b1"]=signal_data["B1"]
        dc3["b2"]=signal_data["B2"]

        l1=file.create_group("lepton_features")
        l1["eta"]=signal_data["lepton_eta1"]
        l1["mass"]=signal_data["lepton_m1"]
        l1["phi"]=signal_data["lepton_phi1"]
        l1["pt"]=signal_data["lepton_pt1"]
        l1["charge"]=signal_data["lepton_charge1"]
        l1["pid"]=signal_data["lepton_pid1"]

        m1=file.create_group("met_features")
        m1["MET"]=signal_data["met_met1"]
        m1["phi"]=signal_data["met_phi1"]
        m1["sumet"]=signal_data["sumet1"]

        j1=file.create_group("jet_features")
        j1["phi"]=jet_data["jet_phi1"]
        j1["pt"]=jet_data["jet_pt1"]
        j1["eta"]=jet_data["jet_eta1"]
        j1["barcode"]=jet_data["jet_barcode1"]
        j1["mass"]=jet_data["jet_m1"]
        j1["btag"]=jet_data["jet_btag1"]

        p1=file.create_group("parton_features")
        p1["phi"]=parton_data["parton_phi1"]
        p1["pt"]=parton_data["parton_pt1"]
        p1["pid"]=parton_data["parton_pid1"]
        p1["eta"]=parton_data["parton_eta1"]
        p1["mass"]=parton_data["parton_m1"]
        p1["barcode"]=parton_data["parton_barcode1"]

        file.create_dataset("events_BDT_feature", data=signal_dataA, dtype=DATA_STRUCTUREA)


    print('Wrote structured data to output file "output_signal.h5".')



if __name__ == "__main__":
    main()
