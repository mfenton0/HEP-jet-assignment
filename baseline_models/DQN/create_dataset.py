from argparse import ArgumentParser
from itertools import permutations

import h5py
import numpy as np
from tqdm import tqdm

import numba


def load_data(input_file: str):
    print(f"Loading HDF5 Data from {input_file}.")

    # Load in data from the concatentated HDF5 file
    with h5py.File(input_file, 'r') as file:
        num_of_jets = file["jet_features"]["num_of_jets"][:]
        MAX_JETS_IN_DATA = num_of_jets.max()

        dataset = np.stack([
            file["jet_features"]["mass"][:, :MAX_JETS_IN_DATA],
            file["jet_features"]["pt"][:, :MAX_JETS_IN_DATA],
            file["jet_features"]["eta"][:, :MAX_JETS_IN_DATA],
            file["jet_features"]["phi"][:, :MAX_JETS_IN_DATA],
            file["jet_features"]["btag"][:, :MAX_JETS_IN_DATA]
        ], axis=-1)

        barcodes = file["jet_features"]["barcode"][:]

        mask = np.zeros((num_of_jets.shape[0], MAX_JETS_IN_DATA), dtype=bool)
        for i, n in enumerate(num_of_jets):
            mask[i, :n] = True

    # Normalize the dataset before hand, this is just for simplicity to avoid doing it later
    # Technically cheating since it means we use testing data statistics as well
    # But we have so many events it doesnt really matter :?
    dataset_mean = dataset[mask].mean(0).reshape(1, 1, -1)
    dataset_std = dataset[mask].std(0).reshape(1, 1, -1)
    dataset = (dataset - dataset_mean) / dataset_std

    return dataset, barcodes, num_of_jets


def create_permutation_dataset(dataset, barcodes, num_of_jets, min_jets, max_jets):
    print(f"Generating Permutation Dataset.")

    # DEFINE TARGET as the following permutation of jets
    target_barcode = np.array([[17, 34, 40, 40]])

    permutation_jets = []
    permutation_barcodes = []
    permutation_targets = []
    event_indices = []

    # Mask out those events that have the desired number of jets.
    jet_mask = (num_of_jets <= max_jets) & (num_of_jets >= min_jets)

    for event_index, (event, event_barcodes, njet) in enumerate(zip(
            tqdm(dataset[jet_mask]),
            barcodes[jet_mask],
            num_of_jets[jet_mask]
    )):
        event = event[:njet]
        event_barcodes = event_barcodes[:njet]

        # Get all possible 4-permutations for this event
        IDX = np.array(list(permutations(range(njet), r=4)))

        # Gather the permutation jet vectors and normalize
        jets = event[IDX]

        # Get each permutations barcode vector and set the target to
        # true if it matches to desired permutation
        jet_barcodes = event_barcodes[IDX]
        target = (jet_barcodes == target_barcode).all(1)

        # Ignore non-fully reconstructable events
        if target.sum() < 2:
            continue

        permutation_jets.append(jets.astype(np.float32))
        permutation_barcodes.append(jet_barcodes)
        permutation_targets.append(target)
        event_indices.extend([event_index] * len(IDX))

    permutation_jets = np.concatenate(permutation_jets)
    permutation_barcodes = np.concatenate(permutation_barcodes)
    permutation_targets = np.concatenate(permutation_targets)
    event_indices = np.array(event_indices)

    return permutation_jets, permutation_barcodes, permutation_targets, event_indices


FUNCTION_TYPE = numba.void(
    numba.float32[:, :, ::1],
    numba.int64[:, ::1],
    numba.int64[::1],
    numba.int64[::1],
    numba.types.ListType(numba.int64[:, ::1]),

    numba.float32[:, :, ::1],
    numba.int64[:, ::1],
    numba.boolean[::1],
    numba.int64[::1],
)


@numba.jit(FUNCTION_TYPE, nopython=True, parallel=True)
def create_permutation_dataset_loop(
        # Inputs
        # ------------------------
        dataset,
        barcodes,
        num_of_jets,
        num_permutations,
        permutation_indices,

        # Outputs
        # ------------------------
        permutation_jets,
        permutation_barcodes,
        permutation_targets,
        event_indices
):
    target_barcode = np.array([[17, 34, 40, 40]])
    permutation_offsets = num_permutations.cumsum() - num_permutations

    for event_index in numba.prange(dataset.shape[0]):
        num_jets = num_of_jets[event_index]
        event_data = dataset[event_index]
        event_barcodes = barcodes[event_index]
        event_offset = permutation_offsets[event_index]

        # Get all possible 4-permutations for this event
        for permutation_index in range(num_permutations[event_index]):
            permutation = permutation_indices[num_jets][permutation_index]

            # Gather the permutation jet vectors and normalize
            current_jets = event_data[permutation]

            # Get each permutations barcode vector and set the target to
            # true if it matches to desired permutation
            current_barcodes = event_barcodes[permutation]
            current_target = (current_barcodes == target_barcode).all()

            output_index = event_offset + permutation_index

            permutation_jets[output_index, :, :] = current_jets
            permutation_barcodes[output_index, :] = current_barcodes
            permutation_targets[output_index] = current_target
            event_indices[output_index] = event_index


def create_permutation_dataset_2(dataset, barcodes, num_of_jets, min_jets, max_jets):
    print(f"Generating Permutation Dataset.")

    # Limit the event selection to full events with the correct number of jets.
    full_event_mask = (
            ((barcodes == 17).sum(1) == 1) &
            ((barcodes == 34).sum(1) == 1) &
            ((barcodes == 40).sum(1) == 2)
    )

    jet_mask = (num_of_jets <= max_jets) & (num_of_jets >= min_jets)
    combined_mask = full_event_mask & jet_mask

    # Select from dataset and convert to final datatype.
    dataset = dataset[combined_mask].astype(np.float32)
    barcodes = barcodes[combined_mask].astype(np.int64)
    num_of_jets = num_of_jets[combined_mask].astype(np.int64)

    # Construct a pre-calculated list of all 4-permutations for each event size.
    permutation_indices = [
        np.atleast_2d(np.array(list(permutations(range(njet), r=4)), dtype=np.int64))
        for njet in range(num_of_jets.max() + 1)
    ]

    # Count the number of permutations for each event.
    num_permutations = np.array([x.shape[0] for x in permutation_indices])
    num_permutations = num_permutations[num_of_jets]

    # Create empty output arrays for the numba function to dump results into.
    num_output_permutations = num_permutations.sum()
    permutation_jets = np.empty((num_output_permutations, 4, 5), dtype=np.float32)
    permutation_barcodes = np.empty((num_output_permutations, 4), dtype=np.int64)
    permutation_targets = np.empty((num_output_permutations,), dtype=bool)
    event_indices = np.empty((num_output_permutations,), dtype=np.int64)

    # Create the dataset inplace
    create_permutation_dataset_loop(
        dataset,
        barcodes,
        num_of_jets,
        num_permutations,
        numba.typed.typedlist.List(permutation_indices),

        permutation_jets,
        permutation_barcodes,
        permutation_targets,
        event_indices
    )

    return permutation_jets, permutation_barcodes, permutation_targets, event_indices


def main(input_file: str, output_file: str, min_jets: int, max_jets: int):
    event_dataset, event_barcodes, num_of_jets = load_data(input_file)

    jets, barcodes, targets, event_indices = create_permutation_dataset_2(
        event_dataset,
        event_barcodes,
        num_of_jets,
        min_jets,
        max_jets
    )

    print(f"Saving result to {output_file}.")
    with h5py.File(output_file, 'w') as file:
        file.create_dataset("jets", data=jets, chunks=(1_000_000, 4, jets.shape[-1]), compression="gzip")
        file.create_dataset("barcodes", data=barcodes, chunks=(1_000_000, 4), compression="gzip")
        file.create_dataset("targets", data=targets, chunks=(1_000_000,), compression="gzip")
        file.create_dataset("event_index", data=event_indices, chunks=(1_000_000,), compression="gzip")


if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument("input_file", type=str)
    parser.add_argument("output_file", type=str)

    parser.add_argument("--min_jets", type=int, default=0)
    parser.add_argument("--max_jets", type=int, default=1_000_000)

    main(**parser.parse_args().__dict__)
