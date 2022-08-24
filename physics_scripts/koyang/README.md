This is the folder to store top mass analysis script before merging in.

1. ttbar_testing_mass_variation.ipynb -- Using SPANET qqb to reconstruct top mass (currently not use)

2. ttbar_testing_mass_variation2.ipynb -- Directly use SPANET prediction of reconstructed top mass to analysis.

3. example-top-ljets.cxx is a copy from https://github.com/mfenton0/HEP-jet-assignment/blob/v4/baseline_models/KLFitter/example-ttbar-ljets.cxx. The difference is I directly read the number of jets and sumet from the file. (Original is set by NOjet and the sum of event.Jet_pt)