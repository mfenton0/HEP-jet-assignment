from glob import glob
from typing import List, Optional
from argparse import ArgumentParser

import h5py
import numpy as np


def extract(file):
    if isinstance(file, h5py.Dataset):
        return file[:]

    database = {}
    for key in file:
        database[key] = extract(file[key])

    return database


def concatenate(head, *tail, path: Optional[List[str]] = None):
    if path is None:
        path = []

    if not isinstance(head, dict):
        return np.concatenate((head, *tail))

    database = {}
    for key in head:
        new_path = path + [key]
        print(f"Concatenating: {'/'.join(new_path)}")
        database[key] = concatenate(head[key], *[d[key] for d in tail], path=new_path)

    return database


def write(input_file, output_file, path: Optional[List[str]] = None):
    if path is None:
        path = []

    for key, value in input_file.items():
        current_subpath = path + [key]
        if isinstance(value, np.ndarray):
            print(f"Creating {'/'.join(current_subpath)}: Shape {value.shape}")
            output_file.create_dataset("/".join(current_subpath), data=value)
        else:
            write(input_file[key], output_file, current_subpath)


def main(input_folder: str, output_file: str):
    print("=" * 40)
    print(f"Reading in files from {input_folder}")
    print("-" * 40)
    global_file = []
    for filename in glob(f"{input_folder}/*.h5"):
        print(f"Reading: {filename}")
        with h5py.File(filename, 'r') as file:
            global_file.append(extract(file))
    print("=" * 40)
    print()

    print("=" * 40)
    print("Concatenating Files")
    print("-" * 40)
    global_file = concatenate(*global_file)
    print("=" * 40)
    print()

    print("=" * 40)
    print(f"Writing Output to {output_file}")
    print("-" * 40)
    with h5py.File(output_file, 'w') as output_file:
        write(global_file, output_file)
    print("=" * 40)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input_folder", type=str, help="Folder of HDF5 files to concatenate.")
    parser.add_argument("output_file", type=str, help="Complete HDF5 file to create for output.")

    args = parser.parse_args()
    main(args.input_folder, args.output_file)
