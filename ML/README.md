# ONNXRuntime
ONNXRuntime enables us to use ML models in python in C++. The following instruction is an example of it with LightGBM. The hyperparameters and training variables in the training phase are defined in the training_configure.yaml file. One can change them by modifying the file. On the other hand, in the inference phase, the feature variables are hard-coded in FwdMatchFunc.C. So, one should check the variable consistency. In near future, an external file, for example, a YAML file, should be used to avoid inconsistency between the training and inference phase.

# How to use
### load O2 and ONNXRuntime
    alienv load O2/latest-dev-o2 ONNXRuntime/latest

### install python modules
    pip install onnxruntime
    pip install onnxmltools
    pip install pyquickhelper
    pip install mlprodict
    pip install skl2onnx
    pip install lightgbm
    pip install docutils
    pip install tf2onnx
    pip install keras2onnx
    pip install awkward

### Simulate a muon per event x 10000 events
    o2-sim -m PIPE MFT MCH MID ABSO SHIL -e TGeant4 -g fwmugen -n 10000 --seed 0

### Reconstruct MFT, MCH and MID
    o2-sim-digitizer-workflow -b
    o2-mft-reco-workflow  --configKeyValues "MFTTracking.forceZeroField=false;MFTTracking.LTFclsRCut=0.0100;" -b
    o2-mch-reco-workflow -b
    o2-mid-digits-reader-workflow -b | o2-mid-reco-workflow -b
    o2-mch-tracks-reader-workflow -b | o2-mid-tracks-reader-workflow -b | o2-muon-tracks-matcher-workflow -b | o2-muon-tracks-writer-workflow -b

### Export training data
    o2-globalfwd-matcher-workflow --configKeyValues "FwdMatching.useMIDMatch=true;FwdMatching.saveMode=2;" -b
    root -q -b FwdMatchTrainingTreeBuilder.C

### Training with LightGBM (export original and onnx format)
    python3 python_training.py training_configure.yaml

### Apply the training 
    o2-globalfwd-matcher-workflow --configKeyValues "FwdMatching.matchFcn=matchExternal;FwdMatching.useMIDMatch=true" -b