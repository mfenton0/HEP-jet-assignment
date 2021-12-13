Implements the permutation-based DNN model described here. 
https://arxiv.org/pdf/1907.11181.pdf

#Requirements
- [torch](https://pytorch.org/get-started/locally/)
- [pytorch_lightning](https://www.pytorchlightning.ai/)
- [numba](https://numba.pydata.org/)
- [h5py](https://pypi.org/project/h5py/)


#Example
``` bash
# Combine datasets created by parsing script 
python concatenate_hdf5.py ./data/raw/semi_leptonic/event_record_top_lep_right_CMS_jetR05_wLeptonFeatures ./data/raw/semi_leptonic/event_record_top_lep_right_CMS_jetR05.h5

# Create the permutation dataset for training DQN
python create_dataset.py ./data/raw/semi_leptonic/event_record_top_lep_right_CMS_jetR05.h5 ./data/baseline/event_record_top_lep_right_CMS_jetR05_4_jets.h5 --min_jets 4

# Run the DQN training (on the gpu, remove --gpu option for cpu training)
# WARNING Very compute intensive
python train.py ./data/baseline/event_record_top_lep_right_CMS_jetR05_4_jets.h5 --logdir baseline_logs --name CMS_4_jets --gpus 1

# Evaluate a trained DQN on the validation dataset (again remove --cuda for cpu evaluation)
python evaluate.py ./baseline_logs/CMS_4_jets/version_0/ --cuda
```
