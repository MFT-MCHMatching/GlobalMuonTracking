// This macro produces a simple ROOT TTree with MFT and MCH track parameters at the matching plane.
// Assumes presence of file `globalfwdparamsatmatchingplane.root` produced by
// `o2-globalfwd-matcher-workflow --configKeyValues "FwdMatching.useMIDMatch=true;FwdMatching.saveMode=2;"`

using MFTTracksType = o2::mft::TrackMFT;
using MCHTracksType = o2::track::TrackParCovFwd;
using MatchesType = o2::dataformats::MatchInfoFwd;
using LabelsType = o2::MCCompLabel;

static const std::size_t sMaxMLFeatures = 50;
float_t mMLInputFeatures[sMaxMLFeatures];
string mMLInputFeaturesName[sMaxMLFeatures];
Int_t Truth;
int mNInputFeatures;

TTree *matchTreeOut = nullptr;

void configOutTree()
{

  mMLInputFeaturesName[0] = "MFT_X";
  mMLInputFeaturesName[1] = "MFT_Y";
  mMLInputFeaturesName[2] = "MFT_Phi";
  mMLInputFeaturesName[3] = "MFT_Tanl";
  mMLInputFeaturesName[4] = "MFT_InvQPt";
  mMLInputFeaturesName[5] = "MFT_Cov00";
  mMLInputFeaturesName[6] = "MFT_Cov01";
  mMLInputFeaturesName[7] = "MFT_Cov11";
  mMLInputFeaturesName[8] = "MFT_Cov02";
  mMLInputFeaturesName[9] = "MFT_Cov12";
  mMLInputFeaturesName[10] = "MFT_Cov22";
  mMLInputFeaturesName[11] = "MFT_Cov03";
  mMLInputFeaturesName[12] = "MFT_Cov13";
  mMLInputFeaturesName[13] = "MFT_Cov23";
  mMLInputFeaturesName[14] = "MFT_Cov33";
  mMLInputFeaturesName[15] = "MFT_Cov04";
  mMLInputFeaturesName[16] = "MFT_Cov14";
  mMLInputFeaturesName[17] = "MFT_Cov24";
  mMLInputFeaturesName[18] = "MFT_Cov34";
  mMLInputFeaturesName[19] = "MFT_Cov44";

  mMLInputFeaturesName[20] = "MCH_X";
  mMLInputFeaturesName[21] = "MCH_Y";
  mMLInputFeaturesName[22] = "MCH_Phi";
  mMLInputFeaturesName[23] = "MCH_Tanl";
  mMLInputFeaturesName[24] = "MCH_InvQPt";
  mMLInputFeaturesName[25] = "MCH_Cov00";
  mMLInputFeaturesName[26] = "MCH_Cov01";
  mMLInputFeaturesName[27] = "MCH_Cov11";
  mMLInputFeaturesName[28] = "MCH_Cov02";
  mMLInputFeaturesName[29] = "MCH_Cov12";
  mMLInputFeaturesName[30] = "MCH_Cov22";
  mMLInputFeaturesName[31] = "MCH_Cov03";
  mMLInputFeaturesName[32] = "MCH_Cov13";
  mMLInputFeaturesName[33] = "MCH_Cov23";
  mMLInputFeaturesName[34] = "MCH_Cov33";
  mMLInputFeaturesName[35] = "MCH_Cov04";
  mMLInputFeaturesName[36] = "MCH_Cov14";
  mMLInputFeaturesName[37] = "MCH_Cov24";
  mMLInputFeaturesName[38] = "MCH_Cov34";
  mMLInputFeaturesName[39] = "MCH_Cov44";

  mMLInputFeaturesName[40] = "MFT_TrackChi2";
  mMLInputFeaturesName[41] = "MFT_NClust";

  mMLInputFeaturesName[42] = "MatchingChi2";

  mNInputFeatures = 43;

  // Config tree
  matchTreeOut = new TTree("matchTree", "MatchTree");
  for (std::size_t nFeature = 0; nFeature < mNInputFeatures; nFeature++)
  {
    matchTreeOut->Branch(mMLInputFeaturesName[nFeature].c_str(), &mMLInputFeatures[nFeature], Form("%s/F", mMLInputFeaturesName[nFeature].c_str()));
  }
  matchTreeOut->Branch("Truth", &Truth, "Truth/I");
}

void loadFeatures(o2::mft::TrackMFT mftTrack, o2::track::TrackParCovFwd mchTrack, o2::dataformats::MatchInfoFwd matchInfo, o2::MCCompLabel label)
{

  mMLInputFeatures[0] = mftTrack.getX();
  mMLInputFeatures[1] = mftTrack.getY();
  mMLInputFeatures[2] = mftTrack.getPhi(),
  mMLInputFeatures[3] = mftTrack.getTanl();
  mMLInputFeatures[4] = mftTrack.getInvQPt();
  mMLInputFeatures[5] = mftTrack.getCovariances()(0, 0);
  mMLInputFeatures[6] = mftTrack.getCovariances()(0, 1);
  mMLInputFeatures[7] = mftTrack.getCovariances()(1, 1);
  mMLInputFeatures[8] = mftTrack.getCovariances()(0, 2);
  mMLInputFeatures[9] = mftTrack.getCovariances()(1, 2);
  mMLInputFeatures[10] = mftTrack.getCovariances()(2, 2);
  mMLInputFeatures[11] = mftTrack.getCovariances()(0, 3);
  mMLInputFeatures[12] = mftTrack.getCovariances()(1, 3);
  mMLInputFeatures[13] = mftTrack.getCovariances()(2, 3);
  mMLInputFeatures[14] = mftTrack.getCovariances()(3, 3);
  mMLInputFeatures[15] = mftTrack.getCovariances()(0, 4);
  mMLInputFeatures[16] = mftTrack.getCovariances()(1, 4);
  mMLInputFeatures[17] = mftTrack.getCovariances()(2, 4);
  mMLInputFeatures[18] = mftTrack.getCovariances()(3, 4);
  mMLInputFeatures[19] = mftTrack.getCovariances()(4, 4);

  mMLInputFeatures[20] = mchTrack.getX();
  mMLInputFeatures[21] = mchTrack.getY();
  mMLInputFeatures[22] = mchTrack.getPhi(),
  mMLInputFeatures[23] = mchTrack.getTanl();
  mMLInputFeatures[24] = mchTrack.getInvQPt();
  mMLInputFeatures[25] = mchTrack.getCovariances()(0, 0);
  mMLInputFeatures[26] = mchTrack.getCovariances()(0, 1);
  mMLInputFeatures[27] = mchTrack.getCovariances()(1, 1);
  mMLInputFeatures[28] = mchTrack.getCovariances()(0, 2);
  mMLInputFeatures[29] = mchTrack.getCovariances()(1, 2);
  mMLInputFeatures[30] = mchTrack.getCovariances()(2, 2);
  mMLInputFeatures[31] = mchTrack.getCovariances()(0, 3);
  mMLInputFeatures[32] = mchTrack.getCovariances()(1, 3);
  mMLInputFeatures[33] = mchTrack.getCovariances()(2, 3);
  mMLInputFeatures[34] = mchTrack.getCovariances()(3, 3);
  mMLInputFeatures[35] = mchTrack.getCovariances()(0, 4);
  mMLInputFeatures[36] = mchTrack.getCovariances()(1, 4);
  mMLInputFeatures[37] = mchTrack.getCovariances()(2, 4);
  mMLInputFeatures[38] = mchTrack.getCovariances()(3, 4);
  mMLInputFeatures[39] = mchTrack.getCovariances()(4, 4);

  mMLInputFeatures[40] = mftTrack.getTrackChi2();
  mMLInputFeatures[41] = mftTrack.getNumberOfPoints();

  mMLInputFeatures[42] = matchInfo.getMFTMCHMatchingChi2();

  Truth = (int)(label.isCorrect());

  return;
}

void FwdMatchTrainingTreeBuilder()
{

  TFile *trainingDataFileIn = new TFile("globalfwdparamsatmatchingplane.root");
  TTree *matchingPlaneParamTree = (TTree *)trainingDataFileIn->Get("MatchingPlaneParams");

  std::vector<MFTTracksType> *mftTracksVec = nullptr;
  matchingPlaneParamTree->SetBranchAddress("mfttracks", &mftTracksVec);

  std::vector<MCHTracksType> *mchTracksVec = nullptr;
  matchingPlaneParamTree->SetBranchAddress("mchtracks", &mchTracksVec);

  std::vector<MatchesType> *matchInfoVec = nullptr;
  matchingPlaneParamTree->SetBranchAddress("matchinfo", &matchInfoVec);

  std::vector<LabelsType> *mcLabelsVec = nullptr;
  matchingPlaneParamTree->SetBranchAddress("MCTruth", &mcLabelsVec);

  auto nTreeEntries = matchingPlaneParamTree->GetEntries();

  if (!nTreeEntries)
  {
    std::cout << " Nothing to process: nTreeEntries = " << nTreeEntries << std::endl;
  }
  std::cout << " Processing " << nTreeEntries << " entries\n";

  TFile *trainingDataFileOut = new TFile("FwdMatchTrainingTree.root", "RECREATE");

  configOutTree();

  for (auto entry = 0; entry < nTreeEntries; entry++)
  {
    std::cout << " Entry # " << entry << std::endl;
    matchingPlaneParamTree->GetEntry(entry);
    auto nMFTTracks = mftTracksVec->size();
    auto nMCHTracks = mchTracksVec->size();
    auto nMatches = matchInfoVec->size();
    auto nLabels = mcLabelsVec->size();

    std::cout << " nMFTTracks = " << nMFTTracks << " ; nMCHTracks = " << nMCHTracks << " ; nMatches = " << nMatches << " ; nLabels " << nLabels << std::endl;

    for (auto nPair = 0; nPair < nMatches; nPair++)
    {
      const auto &mftTrack = mftTracksVec->at(nPair);
      const auto &mchTrack = mchTracksVec->at(nPair);
      const auto &matchInfo = matchInfoVec->at(nPair);
      const auto &label = mcLabelsVec->at(nPair);

      loadFeatures(mftTrack, mchTrack, matchInfo, label);
      matchTreeOut->Fill();
    }
  }
  trainingDataFileOut->Write();
  trainingDataFileOut->Close();
}
