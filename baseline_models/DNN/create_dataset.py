from argparse import ArgumentParser
from itertools import permutations

import h5py
import numpy as np
from tqdm import tqdm

import numba

def load_jet(dataset, max_jets):
    return dataset[:, :max_jets].astype(np.float32)


@numba.jit("void(int64[:, ::1], int64[:, ::1], int64[::1])", nopython=True)
def fill_barcodes(barcodes, targets, target_values):
    for i in range(targets.shape[0]):
        for target_idx in range(targets.shape[1]):
            target = targets[i, target_idx]
            if target >= 0:
                barcodes[i, target] = target_values[target_idx]


@numba.jit(nopython=True)
def compact_masked_data(dataset, barcodes, mask):
    num_jets = mask.sum(1)
    max_real_jets = num_jets.max()
    num_events = dataset.shape[0]
    num_features = dataset.shape[2]

    output_features = np.log(np.zeros((num_events, max_real_jets, num_features), dtype=np.float32))
    output_barcodes = np.zeros((num_events, max_real_jets), dtype=np.int64) - 1

    for i in range(num_events):
        output_features[i, :num_jets[i]] = dataset[i][mask[i]]
        output_barcodes[i, :num_jets[i]] = barcodes[i][mask[i]]

    return output_features, output_barcodes, num_jets


def load_data(input_file: str):
    print(f"Loading HDF5 Data from {input_file}.")

    # Load in data from the concatentated HDF5 file
    with h5py.File(input_file, 'r') as file:
        # Load in the raw data
        pt = file["jet_features"]["pt"][:].astype(np.float32)
        eta = file["jet_features"]["eta"][:].astype(np.float32)
        phi = file["jet_features"]["phi"][:].astype(np.float32)
        mass = file["jet_features"]["mass"][:].astype(np.float32)
        btag = file["jet_features"]["btag"][:].astype(np.float32)

        # Explicit Mask Computation
        mask = (pt > 25) & (np.abs(eta) < 2.5)

        # Construct feature matrix with log transforms on the larger features.
        features = np.stack([np.log(mass + 1), np.log(pt + 1), eta, phi, btag], axis=-1)

        # Construct target barcodes.
        barcodes = np.zeros((features.shape[0], features.shape[1]), dtype=np.int64)

        if "q1" in file["target"]["left_target"].keys():
            targets = np.stack((
                file["target"]["right_target"]["b"][:],
                file["target"]["left_target"]["b"][:],
                file["target"]["left_target"]["q1"][:],
                file["target"]["left_target"]["q2"][:]
            ), axis=1)
        else:
            targets = np.stack((
                file["target"]["left_target"]["b"][:],
                file["target"]["right_target"]["b"][:],
                file["target"]["right_target"]["q1"][:],
                file["target"]["right_target"]["q2"][:]
            ), axis=1)

        target_values = np.array([1, 2, 3, 3])
        fill_barcodes(barcodes, targets, target_values)

    # Normalize the dataset before hand, this is just for simplicity to avoid doing it later
    # Technically cheating since it means we use testing data statistics as well
    # But we have so many events it doesnt really matter :?
    dataset_mean = features[mask].mean(0).reshape(1, 1, -1)
    dataset_std = features[mask].std(0).reshape(1, 1, -1)

    # Disable normalization for btag
    dataset_mean[0, 0, -1] = 0
    dataset_std[0, 0, -1] = 1

    features = (features - dataset_mean) / dataset_std

    # The explicit mask might have holes in the middle, so get rid of them to make dataset compact.
    # This will make computing the permutations much easier later.
    features, barcodes, num_jets = compact_masked_data(features, barcodes, mask)

    return features, barcodes, num_jets


FUNCTION_TYPE = numba.void(
    numba.float32[:, :, ::1],
    numba.int64[:, ::1],
    numba.int64[::1],
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
        event_type,
        permutation_indices,

        # Outputs
        # ------------------------
        permutation_jets,
        permutation_barcodes,
        permutation_targets,
        event_indices
):
    target_barcodes = np.array([
        [1, 2, 3, 3],
        [17, 34, 40, 40],
        [34, 17, 20, 20]
    ])

    permutation_offsets = num_permutations.cumsum() - num_permutations

    for event_index in numba.prange(dataset.shape[0]):
        num_jets = num_of_jets[event_index]
        event_data = dataset[event_index]
        event_barcodes = barcodes[event_index]
        event_offset = permutation_offsets[event_index]

        target_barcode = target_barcodes[event_type[event_index]:event_type[event_index] + 1]

        # Get all possible 4-permutations for this event
        for permutation_index in range(num_permutations[event_index]):
            permutation = permutation_indices[num_jets][permutation_index]

            # Gather the permutation jet vectors and normalize
            current_jets = event_data[permutation]

            # Get each permutations barcode vector and set the target to
            # true if it matches to desired permutation
            current_barcodes = event_barcodes[permutation]
            current_target = (current_barcodes[:4] == target_barcode).all()

            output_index = event_offset + permutation_index

            permutation_jets[output_index, :, :] = current_jets
            permutation_barcodes[output_index, :] = current_barcodes
            permutation_targets[output_index] = current_target
            event_indices[output_index] = event_index


def create_permutation_dataset_2(
        dataset,
        barcodes,
        num_of_jets,
        min_jets,
        max_jets,
        percent_data,
        context_jets: int = 4
):
    print(f"Generating Permutation Dataset.")

    # Limit the event selection to full events with the correct number of jets.
    # right_full_event_mask = (
    #         ((barcodes == 17).sum(1) == 1) &
    #         ((barcodes == 34).sum(1) == 1) &
    #         ((barcodes == 40).sum(1) == 2)
    # )
    #
    # left_full_event_mask = (
    #         ((barcodes == 17).sum(1) == 1) &
    #         ((barcodes == 34).sum(1) == 1) &
    #         ((barcodes == 20).sum(1) == 2)
    # )
    #
    # print(f"Num full left events: {left_full_event_mask.sum()}")
    # print(f"Num full right events: {right_full_event_mask.sum()}")
    #
    # # Mark if events are right or left, important for target selection.
    # event_type = np.zeros(dataset.shape[0], dtype=np.int64)
    # event_type[right_full_event_mask] = 1
    # event_type[left_full_event_mask] = 2
    # assert (event_type > 0).all(), "Events which are not either left or right!"

    # full_event_mask = left_full_event_mask | right_full_event_mask

    event_type = np.zeros(dataset.shape[0], dtype=np.int64)

    full_event_mask = (
            ((barcodes == 1).sum(1) == 1) &
            ((barcodes == 2).sum(1) == 1) &
            ((barcodes == 3).sum(1) == 2)
    )

    jet_mask = (num_of_jets <= max_jets) & (num_of_jets >= min_jets)
    combined_mask = full_event_mask & jet_mask

    if percent_data < 1.0:
        percent_mask = np.random.rand(dataset.shape[0]) < percent_data
        combined_mask &= percent_mask

    print(f"Number of full events: {combined_mask.sum()}")

    # Select from dataset and convert to final datatype.
    dataset = dataset[combined_mask].astype(np.float32)
    barcodes = barcodes[combined_mask].astype(np.int64)
    num_of_jets = num_of_jets[combined_mask].astype(np.int64)
    event_type = event_type[combined_mask]

    # Construct a pre-calculated list of all 4-permutations for each event size.
    permutation_indices = [
        np.atleast_2d(np.array(list(permutations(range(njet), r=context_jets)), dtype=np.int64))
        for njet in range(num_of_jets.max() + 1)
    ]

    # Count the number of permutations for each event.
    num_permutations = np.array([x.shape[0] for x in permutation_indices])
    num_permutations = num_permutations[num_of_jets]

    # Create empty output arrays for the numba function to dump results into.
    num_output_permutations = num_permutations.sum()
    permutation_jets = np.zeros((num_output_permutations, context_jets, 5), dtype=np.float32)
    permutation_barcodes = np.zeros((num_output_permutations, context_jets), dtype=np.int64)
    permutation_targets = np.zeros((num_output_permutations,), dtype=bool)
    event_indices = np.zeros((num_output_permutations,), dtype=np.int64)

    # Create the dataset inplace
    create_permutation_dataset_loop(
        dataset,
        barcodes,
        num_of_jets,
        num_permutations,
        event_type,
        numba.typed.typedlist.List(permutation_indices),

        permutation_jets,
        permutation_barcodes,
        permutation_targets,
        event_indices
    )

    return permutation_jets, permutation_barcodes, permutation_targets, event_indices, combined_mask


def main(input_file: str, output_file: str, min_jets: int, max_jets: int, percent_data: float, context_jets: int):
    event_dataset, event_barcodes, num_of_jets = load_data(input_file)

    jets, barcodes, targets, event_indices, combined_mask = create_permutation_dataset_2(
        event_dataset,
        event_barcodes,
        num_of_jets,
        min_jets,
        max_jets,
        percent_data,
        context_jets
    )

    with h5py.File(input_file, 'r') as file:
        leptons = np.stack([
            np.log(file["lepton_features"]["pt"][:].astype(np.float32) + 1),
            file["lepton_features"]["eta"][:].astype(np.float32),
            file["lepton_features"]["phi"][:].astype(np.float32),
            # np.log(file["met_features"]["sumet"][:, 0].astype(np.float32) + 1),
            np.log(file["met_features"]["MET"][:, 0].astype(np.float32) + 1),
            file["met_features"]["eta"][:, 0].astype(np.float32)
        ], axis=-1)

        leptons = leptons[combined_mask].astype(np.float32)
        leptons -= leptons.mean(0, keepdims=True)
        leptons /= leptons.std(0, keepdims=True)

        leptons = leptons[event_indices]

    CHUNK = 16_000_000

    print(f"Saving result to {output_file}.")
    with h5py.File(output_file, 'w') as file:
        file.create_dataset("jets", data=jets, chunks=(CHUNK, context_jets, jets.shape[-1]), compression="gzip")
        file.create_dataset("leptons", data=leptons, chunks=(CHUNK, leptons.shape[-1]), compression="gzip")
        file.create_dataset("barcodes", data=barcodes, chunks=(CHUNK, context_jets), compression="gzip")
        file.create_dataset("targets", data=targets, chunks=(CHUNK,), compression="gzip")
        file.create_dataset("event_index", data=event_indices, chunks=(CHUNK,), compression="gzip")


if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument("input_file", type=str)
    parser.add_argument("output_file", type=str)

    parser.add_argument("--min_jets", type=int, default=0)
    parser.add_argument("--max_jets", type=int, default=1_000_000)
    parser.add_argument("--context_jets", type=int, default=4)
    parser.add_argument("--percent_data", type=float, default=1.0)

    main(**parser.parse_args().__dict__)
