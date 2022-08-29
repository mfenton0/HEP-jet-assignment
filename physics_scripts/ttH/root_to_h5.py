#!/usr/bin/env python3

import awkward as ak
import h5py
import numpy as np
import tqdm
import uproot


DATA_STRUCTURE = np.dtype(
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
        ("spanet_signal_target", "f4"),
    ]
)


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
    flattened_event["spanet_signal_target"] = event["spanet_signal_target"]

    return flattened_event


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
            "spanet_signal_target",
        ],
    )

    n_evt_tot = min(int(n_evt_max), len(tree)) if n_evt_max else len(tree)
    n_evt_sel = 0

    data = np.zeros(n_evt_tot, dtype=DATA_STRUCTURE)

    print(
        f'Reading {n_evt_tot:,} events from concatenated '
        f'tree "{tree_name}" with {len(tree):,} events.'
    )

    for n_evt in tqdm.tqdm(range(n_evt_tot)):
        event = tree[n_evt]
        #if not passes_selection(event):
        #    continue
        data[n_evt_sel] = get_flattened_event(event)
        n_evt_sel += 1

    print(
        f"Selected {n_evt_sel:,} out of {n_evt_tot:,} events "
        f"({n_evt_sel / n_evt_tot * 100:.1f}%%)."
    )

    return data[:n_evt_sel]


def main():
    """Reads signal and background files and converts them to HDF5."""
   

    signal_files = f"BDTana8.root"
    signal_data = create_structured_data(signal_files, "nominal", 1e7)

    print(f'Created structured data from input files "{signal_files}".')

    with h5py.File("output_signal.h5", "w") as file:
        file.create_dataset("events", data=signal_data, dtype=DATA_STRUCTURE)

    print('Wrote structured data to output file "output_signal.h5".')



if __name__ == "__main__":
    main()
