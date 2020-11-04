
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
#include "MCHSimulation/GeometryTest.h"
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
#endif

#include "MCHTracking/TrackParam.h"
#include "MFTBase/Constants.h"
#include "MFTTracking/Cluster.h"
#include "MFTTracking/Constants.h"
#include "MFTTracking/IndexTableUtils.h"
#include "include/GlobalMuonTrack.h"
#include "include/TrackExtrap.h"
#include "include/tempMCHTrack.h"

#include <iostream>

using MCHTrack = o2::mch::TrackParam;
using MFTTrack = o2::mft::TrackMFT;
using GlobalMuonTrack = o2::track::GlobalMuonTrack;
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
class MUONMatcher {
public:
  MUONMatcher();
  ~MUONMatcher() = default;
  void Clear();
  void SetMatchingPlaneZ(double z) { mMatchingPlaneZ = z; }
  void SetVerbosity(bool v = true) { mVerbose = v; }
  void LoadAbsorber();

  // Track IO
  void loadMCHTracks();
  void loadDummyMCHTracks();
  void loadMFTTracksOut();
  void saveGlobalMuonTracks();
  std::vector<GlobalMuonTrack> getGlobalMuonTracks() const {
    return mGlobalMuonTracks;
  }
  void printMFTLabels() {
    for (auto i = 0; i < (int)mftTrackLabels.getNElements(); i++) {
      for (auto label : mftTrackLabels.getLabels(i)) {
        std::cout << " Track " << i << " label: ";
        label.print();
        std::cout << std::endl;
      }
    }
  }

  void initGlobalTracks(); // Configure Global Tracks with MCH track parameters
  void
  initDummyGlobalTracks(); // Configure Global Tracks with MFT tracks (Dummy)
  void fitTracks();        // Fit all matched tracks

  // Matching methods
  // Position
  double
  matchMFT_MCH_TracksXY(const GlobalMuonTrack &mchTrack,
                        const MFTTrack &mftTrack); // Compute track matching
  //// Position & Angles
  double matchMFT_MCH_TracksXYPhiTanl(const GlobalMuonTrack &mchTrack,
                                      const MFTTrack &mftTrack);
  //// Position, Angles & Charged Momentum
  double matchMFT_MCH_TracksAllParam(const GlobalMuonTrack &mchTrack,
                                     const MFTTrack &mftTrack);
  void setMatchingFunction(double (MUONMatcher::*func)(const GlobalMuonTrack &,
                                                       const MFTTrack &)) {
    mMatchFunc = func;
    if (func == &MUONMatcher::matchMFT_MCH_TracksXY)
      mMatchingHelper.MatchingFunction = "_matchXY";
    if (func == &MUONMatcher::matchMFT_MCH_TracksXYPhiTanl)
      mMatchingHelper.MatchingFunction = "_matchXYPhiTanl";
    if (func == &MUONMatcher::matchMFT_MCH_TracksAllParam)
      mMatchingHelper.MatchingFunction = "_matchAllParams";
    std::cout << " ** MUONMATCHER: Setting matching function => "
              << mMatchingHelper.MatchingFunction << std::endl;
  }
  void setCustomMatchingFunction(double (*func)(const GlobalMuonTrack &,
                                                const MFTTrack &),
                                 std::string nickname) {
    mCustomMatchFunc = func;
    mMatchingHelper.MatchingFunction = nickname;
    std::cout << " ** MUONMATCHER: Setting custom matching function => "
              << mMatchingHelper.MatchingFunction << std::endl;
  }
  void
  runHeavyMatching(); // Finds best match (no search cut, no event separation)
  void runEventMatching(); // Finds best match event-per-event
  void printMatchingPlaneView(int MCHTrackID = 0);
  void exportNMatchingPlaneViews(int nTracks = -1) {
    loadMCHTracks();
    initGlobalTracks();
    if (nTracks < 0 || nTracks > mGlobalMuonTracks.size()) {
      nTracks = mGlobalMuonTracks.size();
    }
    for (auto MCHTrackID = 0; MCHTrackID < nTracks; MCHTrackID++) {
      printMatchingPlaneView(MCHTrackID);
    }
  };

  // Matching cuts
  bool matchingCut(const GlobalMuonTrack &,
                   const MFTTrack &); // Calls configured cut function
  void setCutFunction(bool (MUONMatcher::*func)(const GlobalMuonTrack &,
                                                const MFTTrack &));
  void setCustomCutFunction(bool (*func)(const GlobalMuonTrack &,
                                         const MFTTrack &)) {
    mCustomCutFunc = func;
  }
  //  Built-in cut functions
  bool matchCutDisabled(const GlobalMuonTrack &, const MFTTrack &);
  bool matchCutDistance(const GlobalMuonTrack &, const MFTTrack &);
  bool matchCutDistanceAndAngles(const GlobalMuonTrack &, const MFTTrack &);
  bool matchCutDistanceSigma(const GlobalMuonTrack &, const MFTTrack &);
  bool matchCut3SigmaXYAngles(const GlobalMuonTrack &, const MFTTrack &);
  void setCutParam(int index, double param) {
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
  MCHtoGlobal(MCHTrack &); // Convert MCH Track to GlobalMuonTrack;
  void finalize();

  // Global Muon Track Methods
  void fitGlobalMuonTrack(GlobalMuonTrack &); // Kalman filter
  bool computeCluster(GlobalMuonTrack &, MFTCluster &);
  double (MUONMatcher::*mMatchFunc)(const GlobalMuonTrack &, const MFTTrack &);
  double (*mCustomMatchFunc)(const GlobalMuonTrack &,
                             const MFTTrack &) = nullptr;
  bool (MUONMatcher::*mCutFunc)(const GlobalMuonTrack &, const MFTTrack &);
  bool (*mCustomCutFunc)(const GlobalMuonTrack &, const MFTTrack &) = nullptr;

  // Data Members
  std::vector<MFTTrack> mMFTTracks;
  std::vector<MCHTrack> mMCHTracks;
  std::vector<MFTTrack>
      mMCHTracksDummy; // Dummy MCH tracks at the MFT coordinate system
  std::vector<GlobalMuonTrack> mGlobalMuonTracks;
  std::vector<MFTCluster> mMFTClusters;
  std::vector<int> mtrackExtClsIDs;
  std::vector<o2::itsmft::ROFRecord> mMFTTracksROFs;
  int mNEvents = 0;
  MatchingHelper mMatchingHelper;

  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mftTrackLabels;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mchTrackLabels;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mGlobalTrackLabels;

  // MCH Track Propagation clasee
  o2::mch::TrackExtrap mMCHTrackExtrap;

  double mField_z;
  const double sLastMFTPlaneZ = -77.5;
  double mMatchingPlaneZ = sLastMFTPlaneZ;
  std::vector<double> mCutParams;
  bool mVerbose = false;
  TGeoManager *mGeoManager;
};

//_________________________________________________________________________________________________
std::string getParamString(o2::track::TrackParCovFwd t) {
  std::string param;
  param = " x = " + std::to_string(t.getX()) +
          " y = " + std::to_string(t.getY()) +
          " phi = " + std::to_string(t.getPhi()) +
          " tanl = " + std::to_string(t.getTanl()) +
          " q*pt = " + std::to_string(1.0 / t.getInvQPt());
  return param;
}

//_________________________________________________________________________________________________
std::string getCovString(o2::track::TrackParCovFwd t) {
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
