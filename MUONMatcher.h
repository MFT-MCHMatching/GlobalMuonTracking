
#ifndef MUON_MATCHING
#define MUON_MATCHING

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsPassive/Absorber.h"
#include "DetectorsPassive/Cave.h"
#include "DetectorsPassive/Shil.h"
#include "Field/MagneticField.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTTracking/IOUtils.h"
#include "MFTTracking/TrackCA.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TGeoManager.h"
#include <TGeoGlobalMagField.h>

#include "Math/SMatrix.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#endif

#include "MCHTracking/TrackParam.h"
#include "MFTBase/Constants.h"
#include "MFTTracking/Cluster.h"
#include "MFTTracking/Constants.h"
#include "include/GlobalMuonTrack.h"
#include "include/TrackExtrap.h"
#include "include/tempMCHTrack.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"

#include <iostream>

using MCHTrack = o2::mch::TrackParam;
using MFTTrack = o2::mft::TrackMFT;
using GlobalMuonTrack = o2::track::GlobalMuonTrack;
using GlobalMuonTrackExt = o2::track::GlobalMuonTrackExt;
using MCLabels = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
using MFTCluster = o2::mft::Cluster;

using SMatrix22 = ROOT::Math::SMatrix<double, 2>;
using SMatrix25 = ROOT::Math::SMatrix<double, 2, 5>;
using SMatrix52 = ROOT::Math::SMatrix<double, 5, 2>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using SVector2 = ROOT::Math::SVector<double, 2>;
using SVector4 = ROOT::Math::SVector<double, 4>;
using SVector5 = ROOT::Math::SVector<double, 5>;

using SMatrix44 = ROOT::Math::SMatrix<double, 4>;
using SMatrix45 = ROOT::Math::SMatrix<double, 4, 5>;
using SMatrix54 = ROOT::Math::SMatrix<double, 5, 4>;
using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
using SMatrix55Sym =
  ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

//_________________________________________________________________________________________________
class MUONMatcher
{
 public:
  MUONMatcher();
  ~MUONMatcher() = default;
  void Clear();
  void SetMatchingPlaneZ(double z) { mMatchingPlaneZ = z; }
  void SetVerbosity(bool v = true) { mVerbose = v; }
  void LoadAbsorber();

  // Track IO
  void loadMCHTracks();
  void loadMFTTracksOut();
  void saveGlobalMuonTracks();

  void printMFTLabels()
  {
    for (auto i = 0; i < (int)mftTrackLabels.getNElements(); i++) {
      for (auto label : mftTrackLabels.getLabels(i)) {
        std::cout << " Track " << i << " label: ";
        label.print();
        std::cout << std::endl;
      }
    }
  }

  void initGlobalTracks(); // Configure Global Tracks with MCH track parameters
  void fitTracks();          // Fit all matched tracks

  // Matching methods
  // Position
  double
    matchMFT_MCH_TracksXY(const GlobalMuonTrack& mchTrack,
                          const MFTTrack& mftTrack); // Compute track matching
  //// Position & Angles
  double matchMFT_MCH_TracksXYPhiTanl(const GlobalMuonTrack& mchTrack,
                                      const MFTTrack& mftTrack);
  //// Position, Angles & Charged Momentum
  double matchMFT_MCH_TracksAllParam(const GlobalMuonTrack& mchTrack,
                                     const MFTTrack& mftTrack);

  //// Matching using trained ML
  double matchTrainedML(const GlobalMuonTrack& mchTrack,
                        const MFTTrack& mftTrack);

  void EvaluateML();

  void setMatchingFunction(double (MUONMatcher::*func)(const GlobalMuonTrack&,
                                                       const MFTTrack&))
  {
    mMatchFunc = func;
    if (func == &MUONMatcher::matchMFT_MCH_TracksXY)
      mMatchingHelper.MatchingFunction = "_matchXY";
    if (func == &MUONMatcher::matchMFT_MCH_TracksXYPhiTanl)
      mMatchingHelper.MatchingFunction = "_matchXYPhiTanl";
    if (func == &MUONMatcher::matchMFT_MCH_TracksAllParam)
      mMatchingHelper.MatchingFunction = "_matchAllParams";
    if (func == &MUONMatcher::matchTrainedML)
      mMatchingHelper.MatchingFunction = "_matchML";

    std::cout << " ** MUONMATCHER: Setting matching function => "
              << mMatchingHelper.MatchingFunction << std::endl;
  }
  void setCustomMatchingFunction(double (*func)(const GlobalMuonTrack&,
                                                const MFTTrack&),
                                 std::string nickname)
  {
    mCustomMatchFunc = func;
    mMatchingHelper.MatchingFunction = nickname;
    std::cout << " ** MUONMATCHER: Setting custom matching function => "
              << mMatchingHelper.MatchingFunction << std::endl;
  }
  void setMatchSaveAll(bool val) { mMatchSaveAll = val; }
  void
    runHeavyMatching();    // Finds best match (no search cut, no event separation)
  void runEventMatching(); // Finds best match event-per-event
  bool printMatchingPlaneView(int event, int MCHTrackID = 0);
  void exportNMatchingPlaneViews(int nTracks = -1)
  {
    if (mMatchSaveAll) {
      std::cout << " ** exportNMatchingPlaneViews disabled for option --matchSaveAll" << std::endl;
      return;
    }
    loadMCHTracks();
    initGlobalTracks();
    auto event = 0;
    if (nTracks < 0 or nTracks > mMatchingHelper.nMCHTracks)
      nTracks = mMatchingHelper.nMCHTracks;
    while (nTracks > 0 or event < mNEvents) {
      for (auto& mchTracks : mSortedGlobalMuonTracks) {
        auto MCHTrackID = 0;
        for (auto& mchTrack : mchTracks) {
          printMatchingPlaneView(event, MCHTrackID);
          nTracks--;
          MCHTrackID++;
        }
        event++;
      }
    }
  };

  void exportTrainingDataRoot(int nMCHTracks = -1);

  // Matching cuts
  void disableChargeMatchCut() { mChargeCutEnabled = false; }
  bool matchingCut(const GlobalMuonTrack&,
                   const MFTTrack&); // Calls configured cut function

  double matchingEval(const GlobalMuonTrack&,
                      const MFTTrack&); // Calls configured cut function
  void setCutFunction(bool (MUONMatcher::*func)(const GlobalMuonTrack&,
                                                const MFTTrack&));
  void setCustomCutFunction(bool (*func)(const GlobalMuonTrack&,
                                         const MFTTrack&))
  {
    mCustomCutFunc = func;
  }
  void configureTMVA(std::string filename, float scorecut)
  {
    mTMVAWeightFileName = filename;
    mMLScoreCut = scorecut;
    mTMVAReader = new TMVA::Reader("!Color:!Silent");

    mTMVAReader->AddVariable("MFT_X", &mMCH_MFT_pair[0]);
    mTMVAReader->AddVariable("MFT_Y", &mMCH_MFT_pair[1]);
    mTMVAReader->AddVariable("MFT_Phi", &mMCH_MFT_pair[2]);
    mTMVAReader->AddVariable("MFT_Tanl", &mMCH_MFT_pair[3]);
    mTMVAReader->AddVariable("MFT_InvQPt", &mMCH_MFT_pair[4]);
    mTMVAReader->AddVariable("MFT_Cov00", &mMCH_MFT_pair[5]);
    mTMVAReader->AddVariable("MFT_Cov01", &mMCH_MFT_pair[6]);
    mTMVAReader->AddVariable("MFT_Cov11", &mMCH_MFT_pair[7]);
    mTMVAReader->AddVariable("MFT_Cov02", &mMCH_MFT_pair[8]);
    mTMVAReader->AddVariable("MFT_Cov12", &mMCH_MFT_pair[9]);
    mTMVAReader->AddVariable("MFT_Cov22", &mMCH_MFT_pair[10]);
    mTMVAReader->AddVariable("MFT_Cov03", &mMCH_MFT_pair[11]);
    mTMVAReader->AddVariable("MFT_Cov13", &mMCH_MFT_pair[12]);
    mTMVAReader->AddVariable("MFT_Cov23", &mMCH_MFT_pair[13]);
    mTMVAReader->AddVariable("MFT_Cov33", &mMCH_MFT_pair[14]);
    mTMVAReader->AddVariable("MFT_Cov04", &mMCH_MFT_pair[15]);
    mTMVAReader->AddVariable("MFT_Cov14", &mMCH_MFT_pair[16]);
    mTMVAReader->AddVariable("MFT_Cov24", &mMCH_MFT_pair[17]);
    mTMVAReader->AddVariable("MFT_Cov34", &mMCH_MFT_pair[18]);
    mTMVAReader->AddVariable("MFT_Cov44", &mMCH_MFT_pair[19]);
    mTMVAReader->AddVariable("MCH_X", &mMCH_MFT_pair[20]);
    mTMVAReader->AddVariable("MCH_Y", &mMCH_MFT_pair[21]);
    mTMVAReader->AddVariable("MCH_Phi", &mMCH_MFT_pair[22]);
    mTMVAReader->AddVariable("MCH_Tanl", &mMCH_MFT_pair[23]);
    mTMVAReader->AddVariable("MCH_InvQPt", &mMCH_MFT_pair[24]);
    mTMVAReader->AddVariable("MCH_Cov00", &mMCH_MFT_pair[25]);
    mTMVAReader->AddVariable("MCH_Cov01", &mMCH_MFT_pair[26]);
    mTMVAReader->AddVariable("MCH_Cov11", &mMCH_MFT_pair[27]);
    mTMVAReader->AddVariable("MCH_Cov02", &mMCH_MFT_pair[28]);
    mTMVAReader->AddVariable("MCH_Cov12", &mMCH_MFT_pair[29]);
    mTMVAReader->AddVariable("MCH_Cov22", &mMCH_MFT_pair[30]);
    mTMVAReader->AddVariable("MCH_Cov03", &mMCH_MFT_pair[31]);
    mTMVAReader->AddVariable("MCH_Cov13", &mMCH_MFT_pair[32]);
    mTMVAReader->AddVariable("MCH_Cov23", &mMCH_MFT_pair[33]);
    mTMVAReader->AddVariable("MCH_Cov33", &mMCH_MFT_pair[34]);
    mTMVAReader->AddVariable("MCH_Cov04", &mMCH_MFT_pair[35]);
    mTMVAReader->AddVariable("MCH_Cov14", &mMCH_MFT_pair[36]);
    mTMVAReader->AddVariable("MCH_Cov24", &mMCH_MFT_pair[37]);
    mTMVAReader->AddVariable("MCH_Cov34", &mMCH_MFT_pair[38]);
    mTMVAReader->AddVariable("MCH_Cov44", &mMCH_MFT_pair[39]);

    mTMVAReader->BookMVA("MUONMatcherML", mTMVAWeightFileName);
  }
  //  Built-in cut functions
  bool matchCutDisabled(const GlobalMuonTrack&, const MFTTrack&);
  bool matchCutDistance(const GlobalMuonTrack&, const MFTTrack&);
  bool matchCutDistanceAndAngles(const GlobalMuonTrack&, const MFTTrack&);
  bool matchCutDistanceSigma(const GlobalMuonTrack&, const MFTTrack&);
  bool matchCut3SigmaXYAngles(const GlobalMuonTrack&, const MFTTrack&);
  void setCutParam(int index, double param)
  {
    if (index > ((int)mCutParams.size() - 1))
      mCutParams.resize(index + 1);
    mCutParams[index] = param;
    std::cout << " ** MUONMATCHER: Setting matching cutParam[" << index
              << "] = " << param << std::endl;
  };

 private:
  // Private IO methods
  void loadMFTClusters();

  // Track methods
  GlobalMuonTrack
    MCHtoGlobal(MCHTrack&); // Convert MCH Track to GlobalMuonTrack;
  void finalize();

  // Global Muon Track Methods
  void fitGlobalMuonTrack(GlobalMuonTrack&); // Kalman filter
  bool computeCluster(GlobalMuonTrack&, MFTCluster&);
  double (MUONMatcher::*mMatchFunc)(const GlobalMuonTrack&, const MFTTrack&);
  double (*mCustomMatchFunc)(const GlobalMuonTrack&,
                             const MFTTrack&) = nullptr;
  bool (MUONMatcher::*mCutFunc)(const GlobalMuonTrack&, const MFTTrack&);
  bool (*mCustomCutFunc)(const GlobalMuonTrack&, const MFTTrack&) = nullptr;

  // Data Members
  std::vector<MFTTrack> mMFTTracks;
  std::vector<std::vector<MFTTrack>> mSortedMFTTracks;

  std::vector<std::vector<MCHTrack>> mSortedMCHTracks;

  std::vector<std::vector<GlobalMuonTrack>> mSortedGlobalMuonTracks;

  std::vector<std::vector<GlobalMuonTrackExt>> mSortedGlobalMuonTracksExt;

  std::vector<MFTCluster> mMFTClusters;
  std::vector<int> mtrackExtClsIDs;
  std::vector<o2::itsmft::ROFRecord> mMFTTracksROFs;
  int mNEvents = 0;
  MatchingHelper mMatchingHelper;

  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mftTrackLabels;
  std::vector<std::vector<int>> mftTrackLabelsIDx;

  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>> mSortedMCHTrackLabels;

  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>> mSortedGlobalTrackLabels;

  // MCH Track Propagation clasee
  o2::mch::TrackExtrap mMCHTrackExtrap;

  double mField_z;
  const double sLastMFTPlaneZ = -77.5;
  double mMatchingPlaneZ = sLastMFTPlaneZ;
  std::vector<double> mCutParams;
  bool mVerbose = false;
  TGeoManager* mGeoManager;
  bool mMatchSaveAll = false;
  bool mChargeCutEnabled = true;

  // TMVA interface
  std::string mTMVAWeightFileName;
  TMVA::Reader* mTMVAReader;
  float_t mMCH_MFT_pair[40];
  float_t mMLScoreCut;
};

//_________________________________________________________________________________________________
std::string getParamString(o2::track::TrackParCovFwd t)
{
  std::string param;
  param = " x = " + std::to_string(t.getX()) +
          " y = " + std::to_string(t.getY()) +
          " phi = " + std::to_string(t.getPhi()) +
          " tanl = " + std::to_string(t.getTanl()) +
          " q*pt = " + std::to_string(1.0 / t.getInvQPt());
  return param;
}

//_________________________________________________________________________________________________
std::string getCovString(o2::track::TrackParCovFwd t)
{
  std::string param;
  param = "Cov: (" + std::to_string(t.getCovariances()(0, 0)) + " ; " +
          std::to_string(t.getCovariances()(1, 1)) + " ; " +
          std::to_string(t.getCovariances()(2, 2)) + " ; " +
          std::to_string(t.getCovariances()(3, 3)) + " ; " +
          std::to_string(t.getCovariances()(4, 4)) + ")";
  return param;
}

#include "MUONMatcher.cxx"

#endif /* MUON_MATCHING */
