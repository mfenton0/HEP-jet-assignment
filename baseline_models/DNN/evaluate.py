from glob import glob
from tqdm import tqdm
from argparse import ArgumentParser

import numpy as np
import torch
torch.multiprocessing.set_sharing_strategy("file_system")

from numba import jit, prange
from train import PermutationNetworkModule, PermutationDataModule


def load_checkpoint(log_dir: str, cuda: bool):
    checkpoint = sorted(glob(f"{log_dir}/checkpoints/*"))[-1]
    print(f"Loading Checkpoint: {checkpoint}.")

    checkpoint = torch.load(checkpoint, map_location="cpu")

    options = checkpoint["hyper_parameters"]
    options["gpus"] = int(cuda)

    data_module = PermutationDataModule(options, 8)
    network_module = PermutationNetworkModule(options, data_module)
    network_module.load_state_dict(checkpoint["state_dict"])

    network_module.eval()
    for parameter in network_module.parameters():
        parameter.requires_grad_(False)

    if cuda:
        network_module = network_module.cuda()

    data_module.setup()

    return network_module, data_module


def evaluate_testing_dataset(network_module, data_module, cuda):
    print(f"Running test data through neural network.")
    barcodes = []
    predictions = []
    event_indices = []

    for batch in tqdm(data_module.val_dataloader()):
        data, leptons, barcode, target, event_index = batch
        data = data.view(data.shape[0], -1)
        leptons = leptons.view(data.shape[0], -1)
        features = torch.cat((leptons, data), dim=-1)

        if cuda:
            features = features.cuda()

        prediction = network_module.network(features).cpu()

        barcodes.append(barcode)
        predictions.append(prediction.reshape(-1))
        event_indices.append(event_index)

    barcodes = torch.cat(barcodes)
    predictions = torch.cat(predictions)
    event_indices = torch.cat(event_indices)

    return barcodes, predictions, event_indices


# def extract_event_predictions(barcodes, predictions, event_indices):
#     print(f"Extracting permutation rankings for each event.")
#
#     # First find all events that are present in the testing dataset
#     unique_events = torch.unique(event_indices)
#
#     predicted_barcodes = []
#
#     # For each event, find all of the associated permutations and select the top ranked permutation.
#     for event_index in tqdm(unique_events):
#         event_mask = event_indices == event_index
#         event_barcode = barcodes[event_mask]
#         event_predictions = predictions[event_mask]
#
#         predicted_barcode = event_barcode[event_predictions.argmax()]
#         predicted_barcodes.append(predicted_barcode)
#
#     return torch.stack(predicted_barcodes)

@jit("int64[:, ::1](int64[:, ::1], float32[::1], int64[::1], int64[::1])", nopython=True, parallel=True)
def extract_event_predictions_loop(barcodes, predictions, event_indices, bounds):
    num_events = bounds.shape[0] - 1
    context_jets = barcodes.shape[1]
    predicted_barcodes = np.zeros((num_events, context_jets), dtype=np.int64)

    for i in prange(num_events):
        lower = bounds[i]
        upper = bounds[i + 1]

        event_barcode = barcodes[lower:upper]
        event_predictions = predictions[lower:upper]
        predicted_barcodes[i] = event_barcode[event_predictions.argmax()]

    return predicted_barcodes


def extract_event_predictions(barcodes, predictions, event_indices):
    print(f"Extracting permutation rankings for each event.")

    event_indices = event_indices.numpy()
    predictions = predictions.numpy()
    barcodes = barcodes.numpy()

    bounds = np.unique(event_indices, return_index=True)[1]
    bounds = np.append(bounds, len(event_indices) + 1)

    return torch.from_numpy(extract_event_predictions_loop(barcodes, predictions, event_indices, bounds))


def main(log_dir: str, cuda: bool):
    network_module, data_module = load_checkpoint(log_dir, cuda)
    barcodes, predictions, event_indices = evaluate_testing_dataset(network_module, data_module, cuda)
    predicted_barcodes = extract_event_predictions(barcodes, predictions, event_indices)

    uniform_target_barcode = torch.tensor([[1, 2, 3, 3]])
    uniform_accuracies = predicted_barcodes[:, :4] == uniform_target_barcode

    left_target_barcode = torch.tensor([[17, 34, 40, 40]])
    left_accuracies = predicted_barcodes[:, :4] == left_target_barcode

    right_target_barcode = torch.tensor([[34, 17, 20, 20]])
    right_accuracies = predicted_barcodes[:, :4] == right_target_barcode

    accuracies = left_accuracies | right_accuracies | uniform_accuracies

    all_accuracy = accuracies.all(1).float().mean().item()
    b_lep_accuracy = accuracies[:, 0].float().mean().item()
    b_had_accuracy = accuracies[:, 1].float().mean().item()
    W_accuracy = accuracies[:, 2:].all(1).float().mean().item()

    print(f"All Accuracy: {100 * all_accuracy:.2f}%.")
    print(f"Leptonic b Quark Accuracy: {100 * b_lep_accuracy:.2f}%.")
    print(f"Hadronic b Quark Accuracy: {100 * b_had_accuracy:.2f}%.")
    print(f"W Boson Accuracy: {100 * W_accuracy:.2f}%.")

    torch.save([event_indices, accuracies], "accuracies.torch")


if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument("log_dir", type=str,
                        help="The output directory during training, typically ..._logs/version_0/")

    parser.add_argument("--cuda", action='store_true')

    main(**parser.parse_args().__dict__)
