HEP jet assignment - Data preparation(main page)
===

## Abstract 
This is a repository for the jet assignment project using state-of-the-art Machine Learning method.

There is two main part in the repository, `madgraph` and `analysis_script`. The `madgraph` folder contains the configuration and auto-run script for generating Monte Carlo simulation data.

## Madgraph
In this project, we generate the data base on the follwing model.

1. Fully hadronic top decay[[link]](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/madgraph/pptt_preparation):
$\qquad p\quad p\quad \to\quad t\quad \bar{t}\quad \to\quad W^{+}b\quad W^{-}\bar{b}\quad \to q_{1}q_{2}bq_{3}q_{4}\bar{b}$
2. Standard Model Higgs boson produced in association with top quarks[[link]](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/madgraph/ppttH_preparation):
$\qquad p\quad p\quad \to\quad t\quad \bar{t}\quad H\quad \to\quad W^{+}b\quad W^{-}\bar{b}\quad b\bar{b}\to q_{1}q_{2}b_{1}q_{3}q_{4}\bar{b}_{1}\quad b_{2}\bar{b}_{2}$
3. Four top production(fully hadronic decay)[[link]](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/madgraph/four_top_preparation):
$\qquad p\quad p\quad \to\quad t\quad \bar{t}\quad t\quad \bar{t}\quad \to\quad W^{+}_{1}b_{1}W^{+}_{2}b_{2} W^{-}_{1}\bar{b}_{1}W^{-}_{2}\bar{b}_{2}\quad \to q_{1}q_{2}b_{1}q_{3}q_{4}\bar{b}_{1}\quad q_{5}q_{6}b_{2}q_{7}q_{8}\bar{b}_{2}$
4. Semi-leptonic top decay[[link]](https://github.com/davidho27941/HEP-jet-assignment/tree/v2/madgraph/ttbar-semi-lep_preparation):
$\qquad p\quad p\quad \to\quad t\quad \bar{t}\quad \to\quad W^{+}b\quad W^{-}\bar{b}\quad \to q_{1}q_{2}bl^{-}\nu\bar{b}$


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

