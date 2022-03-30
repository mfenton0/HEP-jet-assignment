from argparse import ArgumentParser
from typing import List, Optional

import h5py

import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader

import pytorch_lightning as pl
from pytorch_lightning.loggers import TensorBoardLogger

# Options taken from https://arxiv.org/pdf/1907.11181.pdf page 12
OPTIONS_4_JETS = {
    "dataset": "./data/baseline/semi_leptonic_CMS_permutation_4_jets.h5",
    "train_test_split": 0.8,
    "batch_size": 12_000,

    "hidden_dims": [512, 256, 128, 64, 32],
    "learning_rate": 0.01,
    "l2_norm": 1e-8,
    "dropout": 0.0
}

OPTIONS_5_JETS = {
    "dataset": "./data/baseline/semi_leptonic_CMS_permutation_5_jets.h5",
    "train_test_split": 0.8,
    "batch_size": 30_000,

    "hidden_dims": [512, 256, 128, 64, 32, 16],
    "learning_rate": 0.01,
    "l2_norm": 1e-8,
    "dropout": 0.0
}

OPTIONS_6_JETS = {
    "dataset": "./data/baseline/semi_leptonic_CMS_permutation_6_jets.h5",
    "train_test_split": 0.8,
    "batch_size": 30_000,

    "hidden_dims": [512, 256, 128, 64, 32],
    "learning_rate": 0.01,
    "l2_norm": 1e-8,
    "dropout": 0.0
}


class Network(nn.Sequential):
    def __init__(self, input_dim: int, hidden_dims: List[int], dropout: float = 0.0):
        layers = []

        current_dim = input_dim
        for next_dim in hidden_dims:
            layers.append(nn.Linear(current_dim, next_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))

            current_dim = next_dim

        layers.append(nn.Linear(current_dim, 1))

        super(Network, self).__init__(*layers)


# noinspection PyAbstractClass
class PermutationDataModule(pl.LightningDataModule):
    def __init__(self, options, num_workers):
        super().__init__()
        self.options = options
        self.num_workers = num_workers
        self.dataset = options["dataset"]
        self.batch_size = options["batch_size"]
        self.train_test_split = options["train_test_split"]
        self.testing_file = options.get("testing_file", None)

        with h5py.File(self.dataset, 'r') as file:
            self.context_jets, self.jet_features = file["jets"].shape[1:]
            self.lepton_features = file["leptons"].shape[1]
            self.flat_features = self.context_jets * self.jet_features + self.lepton_features

    def setup(self, stage=None):
        with h5py.File(self.dataset, 'r') as file:
            jets = file["jets"][:]
            leptons = file["leptons"][:]
            barcodes = file["barcodes"][:]
            targets = file["targets"][:]
            event_index = file["event_index"][:]

        train_idx = round((event_index.max() + 1) * self.train_test_split)

        train_mask = event_index < train_idx
        permutation_train_idx = train_mask.sum()

        # noinspection PyAttributeOutsideInit
        self.training_dataset = TensorDataset(
            torch.from_numpy(jets[:permutation_train_idx]).float(),
            torch.from_numpy(leptons[:permutation_train_idx]).float(),
            torch.from_numpy(barcodes[:permutation_train_idx]),
            torch.from_numpy(targets[:permutation_train_idx]),
            torch.from_numpy(event_index[:permutation_train_idx])
        )

        # noinspection PyAttributeOutsideInit
        self.validation_dataset = TensorDataset(
            torch.from_numpy(jets[permutation_train_idx:]).float(),
            torch.from_numpy(leptons[permutation_train_idx:]).float(),
            torch.from_numpy(barcodes[permutation_train_idx:]),
            torch.from_numpy(targets[permutation_train_idx:]),
            torch.from_numpy(event_index[permutation_train_idx:])
        )

    def train_dataloader(self):
        return DataLoader(
            self.training_dataset,
            batch_size=self.batch_size,
            shuffle=True,
            pin_memory=self.options["gpus"] > 0,
            num_workers=self.num_workers,
            persistent_workers=True
        )

    def val_dataloader(self):
        return DataLoader(
            self.validation_dataset,
            batch_size=self.batch_size,
            pin_memory=self.options["gpus"] > 0,
            num_workers=self.num_workers,
            persistent_workers=True
        )


# noinspection PyAbstractClass
class PermutationNetworkModule(pl.LightningModule):
    def __init__(self, options, data_module: PermutationDataModule):
        super(PermutationNetworkModule, self).__init__()

        self.options = options
        self.save_hyperparameters(options)

        self.network = Network(
            data_module.flat_features,
            options["hidden_dims"],
            options["dropout"]
        )

        self.loss = nn.BCEWithLogitsLoss()

    def training_step(self, batch, batch_idx):
        data, leptons, barcodes, targets, event_index = batch
        data = data.view(data.shape[0], -1)
        leptons = leptons.view(data.shape[0], -1)
        features = torch.cat((leptons, data), dim=-1)

        predictions = self.network(features)

        loss = self.loss(predictions.squeeze(), targets.float())
        accuracy = ((predictions > 0) == targets).float().mean()

        self.log('train_loss', loss)
        self.log('train_accuracy', accuracy)

        return loss

    def validation_step(self, batch, batch_idx):
        data, leptons, barcodes, targets, event_index = batch
        data = data.view(data.shape[0], -1)
        leptons = leptons.view(data.shape[0], -1)
        features = torch.cat((leptons, data), dim=-1)

        predictions = self.network(features)

        loss = self.loss(predictions.squeeze(), targets.float())
        accuracy = ((predictions > 0) == targets).float().mean()

        self.log('val_loss', loss)
        self.log('val_accuracy', accuracy)

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(
            self.parameters(),
            lr=self.options["learning_rate"],
            weight_decay=self.options["l2_norm"]
        )
        return optimizer


def parse_arguments():
    parser = ArgumentParser()

    parser.add_argument("dataset", type=str)
    parser.add_argument("--batch_size", type=int, default=None)

    parser.add_argument("--gpus", type=int, default=0)
    parser.add_argument("--logdir", type=str, default="baseline_logs")
    parser.add_argument("--name", type=str, default=None)
    parser.add_argument("--epochs", type=int, default=1_000)
    parser.add_argument("--num_workers", type=int, default=8)

    return parser.parse_args()


def main(
        dataset: str,
        logdir: str,
        name: Optional[str],
        gpus: int,
        batch_size: Optional[int],
        epochs: int,
        num_workers: int
):
    if "4_jets" in dataset:
        print("Training on >= 4 jets.")
        options = OPTIONS_4_JETS.copy()
    elif "5_jets" in dataset:
        print("Training on >= 5 jets.")
        options = OPTIONS_5_JETS.copy()
    elif "6_jets" in dataset:
        print("Training on >= 6 jets.")
        options = OPTIONS_6_JETS.copy()
    else:
        options = OPTIONS_4_JETS.copy()

    options["dataset"] = dataset
    options["gpus"] = gpus

    if batch_size is not None:
        options["batch_size"] = batch_size

    data_module = PermutationDataModule(options, num_workers)
    network_module = PermutationNetworkModule(options, data_module)

    logger = TensorBoardLogger(save_dir=logdir, name=name)

    trainer = pl.Trainer(
        logger=logger,
        gpus=gpus,
        weights_summary='full',
        max_epochs=epochs,
        strategy="ddp" if gpus > 1 else None
    )

    trainer.fit(network_module, data_module)


if __name__ == "__main__":
    arguments = parse_arguments()
    main(**arguments.__dict__)
