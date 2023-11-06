* Datasets (JSF = 0p96, 0p98, 1p00, 1p02, 1p04)
    * [ttbar_testing_final_mass_variation_JSF={JSF}](https://cernbox.cern.ch/s/s1LD3rNYr1P878a) -- The folders for dataset with different top mass and JSF for templates (For JSF = 1, the KLFitter applies mt fixed method, needs to use another file. Still looking for the link of the file.)
    * [ttbar_testing_pseudodata_JSF{JSF}_with_spanet_pdnn_KLFitter.h5](https://cernbox.cern.ch/s/koUoDJfV3ojw1J1) -- The datset for final mass extraction

* Data Analysis
    * ttbar_testing_mass_variation_feature.ipynb -- Analysis for cuts and number of jets

* 1D Analysis
Using only the dataset with JSF = 1.00
    * Python Notebook
        * ttbar_testing_mass_variation_top_mass_with_true.ipynb -- Template construction and calibration

* 2D Analysis
    * Python Notebook
        * ttbar_testing_mass_variation_top_mass_JSF_with_true.ipynb -- Template construction and calibration
        * ttbar_testing_mass_variation_plots.ipynb -- Complementary plots for publication
        * ttbar_testing_unknown_mass_top_mass_JSF.ipynb -- Final 2D mass extraction