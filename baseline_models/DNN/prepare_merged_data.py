from glob import glob
from typing import List, Optional
from argparse import ArgumentParser

import h5py
import numpy as np

from concatenate_data import extract, write


def main(input_folder: str, output_folder: str):
    for input_filename in glob(f"{input_folder}/*.h5"):
        print(f"Reading: {input_filename}")
        with h5py.File(input_filename, 'r') as input_file:
            data = extract(input_file)

        # Make all events look like right events, for simplicity of code in the future.
        if "lepton" in data["target"]["left_target"]:
            data["target"]["left_target"], data["target"]["right_target"] = data["target"]["right_target"], data["target"]["left_target"]

        output_filename = input_filename.replace(input_folder, output_folder)
        print(f"Writing: {output_filename}")
        with h5py.File(output_filename, 'w') as output_file:
            write(data, output_file, verbose=False)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input_folder", type=str, help="Folder of HDF5 files to concatenate.")
    parser.add_argument("output_file", type=str, help="Complete HDF5 file to create for output.")

    args = parser.parse_args()
    main(args.input_folder, args.output_file)
