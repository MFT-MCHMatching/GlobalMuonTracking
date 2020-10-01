
#ifndef MUON_MATCHING
#define MUON_MATCHING

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TGeoGlobalMagField.h>
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "MFTTracking/TrackCA.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "MFTBase/GeometryTGeo.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "MFTTracking/IOUtils.h"
#include "DataFormatsITSMFT/ROFRecord.h"

#include "TFile.h"
#include "TTree.h"
#include "Math/SMatrix.h"
#endif

#include "include/TrackExtrap.h"
#include "include/GlobalMuonTrack.h"
#include "include/tempMCHTrack.h"
#include "MCHTracking/TrackParam.h"
#include "MFTTracking/Cluster.h"
#include "MFTTracking/Constants.h"
#include "MFTBase/Constants.h"
#include "MFTTracking/IndexTableUtils.h"



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
using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;




class MUONMatcher
{
public:
  MUONMatcher();
  ~MUONMatcher() = default;
  void Clear();
  void SetMatchingPlaneZ(double z) { mMatchingPlaneZ = z;}
  void SetVerbosity(bool v = true) { mVerbose = v; }


  // Track IO
  void loadMCHTracks();
  void loadDummyMCHTracks();
  void loadMFTTracksOut();
  void saveGlobalMuonTracks();
  std::vector<GlobalMuonTrack> getGlobalMuonTracks() const { return mGlobalMuonTracks;}
  void printMFTLabels() {
    for ( auto i = 0 ; i < mftTrackLabels.getNElements() ; i++) {
      for ( auto label : mftTrackLabels.getLabels(i) )
      {
        std::cout << " Track " << i << " label: "; label.print(); std::cout << std::endl;
      }
    }
  }

  void initGlobalTracks(); // Configure Global Tracks with MCH track parameters
  void initDummyGlobalTracks(); // Configure Global Tracks with MFT tracks (Dummy)
  void fitTracks(); //Fit all matched tracks


  // Matching methods
  // Position
  double matchMFT_MCH_TracksXY(const GlobalMuonTrack& mchTrack, const MFTTrack& mftTrack); // Compute track matching
  //// Position & Angles
  double matchMFT_MCH_TracksXYPhiTanl(const GlobalMuonTrack& mchTrack, const MFTTrack& mftTrack);
  //// Position, Angles & Charged Momentum
  double matchMFT_MCH_TracksFull(const GlobalMuonTrack& mchTrack, const MFTTrack& mftTrack);
  void setMatchingFunction(double (MUONMatcher::*func)(const GlobalMuonTrack&, const MFTTrack&)) { mMatchFunc = func; }
  void setCustomMatchingFunction(double (*func)(const GlobalMuonTrack&, const MFTTrack&)) { mCustomMatchFunc = func; }
  void runHeavyMatching(); // Finds best match (no search cut, no event separation)
  void runEventMatching(); // Finds best match event-per-event

  // Matching cuts
  bool matchingCut(const GlobalMuonTrack&, const MFTTrack&); // Calls configured cut function
  void setCutFunction(bool (MUONMatcher::*func)(const GlobalMuonTrack&, const MFTTrack&)) { mCutFunc = func; }
  void setCustomCutFunction(bool (*func)(const GlobalMuonTrack&, const MFTTrack&)) { mCustomCutFunc = func; }
  //  Built-in cut functions
  bool matchCutDisabled(const GlobalMuonTrack&, const MFTTrack&);
  bool matchCutDistance(const GlobalMuonTrack&, const MFTTrack&);
  bool matchCutDistanceSigma(const GlobalMuonTrack&, const MFTTrack&);
  void setCutDistanceParam(double distance) { mCutDistanceParam = distance; };



private:

  // Private IO methods
  void loadMFTClusters();

  // Track methods
  GlobalMuonTrack MCHtoGlobal(MCHTrack&); // Convert MCH Track to GlobalMuonTrack;
  void ComputeLabels();


  // Global Muon Track Methods
  void fitGlobalMuonTrack(GlobalMuonTrack&); // Kalman filter
  bool computeCluster(GlobalMuonTrack&, MFTCluster&);
  double (MUONMatcher::*mMatchFunc)(const GlobalMuonTrack&, const MFTTrack&) ;
  double (*mCustomMatchFunc)(const GlobalMuonTrack&, const MFTTrack&) = nullptr;
  bool (MUONMatcher::*mCutFunc)(const GlobalMuonTrack&, const MFTTrack&) ;
  bool (*mCustomCutFunc)(const GlobalMuonTrack&,const MFTTrack&) = nullptr;



  // Data Members
  std::vector<MFTTrack> mMFTTracks;
  std::vector<MCHTrack> mMCHTracks;
  std::vector<MFTTrack> mMCHTracksDummy; // Dummy MCH tracks at the MFT coordinate system
  std::vector<GlobalMuonTrack> mGlobalMuonTracks;
  std::vector<MFTCluster> mMFTClusters;
  std::vector<int> mtrackExtClsIDs;
  std::vector<o2::itsmft::ROFRecord> mMFTTracksROFs;
  std::vector<int> mFakeMatches;
  std::vector<int> mGoodMatches;
  int mTotalFakeMatches = 0;
  int mTotalGoodMatches = 0;
  int mNEvents = 0;


  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mftTrackLabels;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mchTrackLabels;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mGlobalTrackLabels;


  // MCH Track Propagation clasee
  o2::mch::TrackExtrap mMCHTrackExtrap;

  double mField_z;
  const double sLastMFTPlaneZ = -77.5;
  double mMatchingPlaneZ = sLastMFTPlaneZ;
  double mCutDistanceParam = 1.0;
  bool mVerbose = false;

};

#include "MUONMatcher.cxx"


#endif /* MUON_MATCHING */
