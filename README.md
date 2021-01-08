HEP jet assignment - Data preparation(main page)
===

## Abstract 
This is a repository for the jet assignment project using state-of-the-art Machine Learning method.

There is two main part in the repository, `madgraph` and `analysis_script`. The `madgraph` folder contains the configuration and auto-run script for generating Monte Carlo simulation data.

## Madgraph
In this project, we generate the data base on the follwing model.

1. Fully hadronic top decay[[link]](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/madgraph/pptt_preparation):  
```p p > t t~ QED=0, (t > W+ b, W+ > j j), (t~ > w- b~, w- > j j )```  
2. Standard Model Higgs boson produced in association with top quarks[[link]](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/madgraph/ppttH_preparation):  
```p p > t t~ h , (t > W+ b, W+ > j j), (t~ > w- b~, w- > j j ), (h > b b~ )```  
3. Four top production(fully hadronic decay)[[link]](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/madgraph/four_top_preparation):  
```p p > t t~ t t~ QED=0, (t > W+ b, W+ > j j), (t~ > w- b~, w- > j j )```  
4. Semi-leptonic top decay[[link]](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/madgraph/ttbar-semi-lep_preparation):  
```p p > t t~ QED=0, (t > W+ b, W+ > j j), (t~ > W- b~, W- > l- vl~)```  

## Analysis 

The script for analysis events can be found in this [folder](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/analysis_script).

The supported analysis method in this repository is:
1. Delta R matching(truth matching)
2. $\chi^{2}$ reconstruction(Only available for two models[^1])
3. Cutflow[^2]
4. Gaussian fitting for finding $\sigma$ for reconstructed invariant mass. 



[^1]: Fully hadronic top decay and Standard Model Higgs boson produced in association with top quarks
[^2]: Only support number of cuts lager than 2 and less than 6.
###### tags: `Particle Physics`, `Machine Learning`, `Top quark`, `Transformer`, `SPA-Net`, `SPAttER`

