// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "onnxruntime/core/session/experimental_onnxruntime_cxx_api.h"
#include "GlobalTracking/MatchGlobalFwd.h"

using o2::dataformats::GlobalFwdTrack;
using o2::globaltracking::CutFunc_t;
using o2::globaltracking::MatchingFunc_t;
using o2::track::TrackParCovFwd;

std::vector<float> getVariables(const GlobalFwdTrack &mchTrack, const TrackParCovFwd &mftTrack)
{
  Float_t MFT_X      = mftTrack.getX();
  Float_t MFT_Y      = mftTrack.getY();
  Float_t MFT_Phi    = mftTrack.getPhi();
  Float_t MFT_Tanl   = mftTrack.getTanl();
  Float_t MFT_InvQPt = mftTrack.getInvQPt();

  Float_t MCH_X      = mchTrack.getX();
  Float_t MCH_Y      = mchTrack.getY();
  Float_t MCH_Phi    = mchTrack.getPhi();
  Float_t MCH_Tanl   = mchTrack.getTanl();
  Float_t MCH_InvQPt = mchTrack.getInvQPt();
  
  Float_t MFT_TrackChi2   = mftTrack.getTrackChi2();
        
  Float_t Delta_X      = MFT_X      - MCH_X;
  Float_t Delta_Y      = MFT_Y      - MCH_Y;
  Float_t Delta_Phi    = MFT_Phi    - MCH_Phi;
  Float_t Delta_Tanl   = MFT_Tanl   - MCH_Tanl;
  Float_t Delta_InvQPt = MFT_InvQPt - MCH_InvQPt;

  std::vector<float> input_tensor_values{
      MFT_X,
      MFT_Y,
      MFT_Phi,
      MFT_Tanl,
      MFT_InvQPt,
      MCH_X,
      MCH_Y,
      MCH_Phi,
      MCH_Tanl,
      MCH_InvQPt,
      MFT_TrackChi2,
      Delta_X,
      Delta_Y,
      Delta_Phi,
      Delta_Tanl,
      Delta_InvQPt
  };

  return input_tensor_values;  
}

std::string model_file = "model.onnx";
  
Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "example-model-explorer");
Ort::SessionOptions session_options;
Ort::Experimental::Session session = Ort::Experimental::Session(env, model_file, session_options);  // access experimental components via the Experimental namespace

MatchingFunc_t matchONNX = [](const GlobalFwdTrack &mchTrack, const TrackParCovFwd &mftTrack) -> double
{  

  std::vector<std::string> input_names;
  std::vector<std::vector<int64_t>> input_shapes;
  std::vector<std::string> output_names;
  std::vector<std::vector<int64_t>> output_shapes;
  
  input_names = session.GetInputNames();
  input_shapes = session.GetInputShapes();
  output_names = session.GetOutputNames();
  output_shapes = session.GetOutputShapes();
  
  auto input_shape = input_shapes[0];
  input_shape[0] = 1;
  
  std::vector<float> input_tensor_values;
  input_tensor_values = getVariables(mchTrack,mftTrack);

  std::vector<Ort::Value> input_tensors;
  input_tensors.push_back(Ort::Experimental::Value::CreateTensor<float>
  			  (input_tensor_values.data(), input_tensor_values.size(), input_shape));
  
  std::vector<Ort::Value> output_tensors = session.Run(input_names, input_tensors, output_names);

  const float* output_value = output_tensors[0].GetTensorData<float>();
  
  auto score = 1-output_value[0];

  return score;
};

MatchingFunc_t *getMatchingFunction()
{
  return &matchONNX;
}

CutFunc_t sameCharge = [](const GlobalFwdTrack &mchTrack, const TrackParCovFwd &mftTrack) -> bool
{
  // Example of external MFTMCH candidate cut function
  return mchTrack.getCharge() == mftTrack.getCharge();
};

CutFunc_t *getCutFunction()
{
  return &sameCharge;
}
