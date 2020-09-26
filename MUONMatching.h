
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

class MUONMatching
{
public:
  MUONMatching();
  ~MUONMatching() = default;
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


  // Matching methods
  // Position
  double matchMFT_MCH_TracksXY(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack); // Compute track matching
  //// Position & Angles
  double matchMFT_MCH_TracksXYPhiTanl(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack);
  //// Position, Angles & Charged Momentum
  double matchMFT_MCH_TracksFull(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack);
  void setMatchingFunction(double (MUONMatching::*func)(GlobalMuonTrack&, MFTTrack&)) { mMatchFunc = func; }
  void setCustomMatchingFunction(double (*func)(GlobalMuonTrack&, MFTTrack&)) { mCustomMatchFunc = func; }
  void runHeavyMatching(); // Finds best match (no search window, no event separation)
  void runEventMatching(); // Finds best match event-per-event

  void fitTracks(); //Fit all matched tracks


private:

  // Private IO methods
  void loadMFTClusters();

  // Track methods
  GlobalMuonTrack MCHtoGlobal(MCHTrack&); // Convert MCH Track to GlobalMuonTrack;
  void ComputeLabels();


  // Global Muon Track Methods
  void fitGlobalMuonTrack(GlobalMuonTrack&); // Kalman filter
  bool computeCluster(GlobalMuonTrack&, MFTCluster&);
  double (MUONMatching::*mMatchFunc)(GlobalMuonTrack&,MFTTrack&) ;
  double (*mCustomMatchFunc)(GlobalMuonTrack&,MFTTrack&) = nullptr;

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
  bool mVerbose = false;

};

#include "MUONMatching.cxx"


#endif /* MUON_MATCHING */
