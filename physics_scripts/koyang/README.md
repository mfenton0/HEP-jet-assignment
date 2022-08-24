This is the folder to store top mass analysis script before merging in. Recommand gdown to download the datasets.

1. ttbar_testing_mass_variation.ipynb -- Using SPANET qqb to reconstruct top mass (currently not use)

2. ttbar_testing_mass_variation2.ipynb -- Directly use SPANET prediction of reconstructed top mass to analysis. Uses top mass = 173 as sudo experiment.
## ttbar_testing_mass_variation: https://drive.google.com/u/0/uc?id=1qMBHbaLKLn3hO_q7WJzoQPYrsI3mMAJu&export=download

3. ttbar_testing_mass_variation2.ipynb -- Different from 2, uses ttbar_testing_sm as sudo experiment.
## ttbar_testing_mass_variation: https://drive.google.com/u/0/uc?id=1qMBHbaLKLn3hO_q7WJzoQPYrsI3mMAJu&export=download
## ttbar_testing_sm: https://drive.google.com/u/0/uc?id=1IwHJxVwMP8it75lQ0ft88bCj26Ci2MdT&export=download

example-top-ljets.cxx is a copy from https://github.com/mfenton0/HEP-jet-assignment/blob/v4/baseline_models/KLFitter/example-ttbar-ljets.cxx. The difference is I directly read the number of jets and sumet from the file. (Original is set by NOjet and the sum of event.Jet_pt)