
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
   void SetMatchingPlane(double z) { mMatchingPlaneZ = z;}

   // Track IO
   void loadMCHTracks();
   void loadDummyMCHTracks();
   void loadMFTTracksOut();
   void saveGlobalMuonTracks();
   std::vector<GlobalMuonTrack> getGlobalMuonTracks() const { return mGlobalMuonTracks;}

   void loadROFrameData(int); // Loads data from a ROFrame
   void initGlobalTracks(); // Configure Global Tracks with MCH track parameters
   void initDummyGlobalTracks(); // Configure Global Tracks with MFT tracks (Dummy)


   // Matching methods
   // Position
   GlobalMuonTrack matchMFT_MCH_TracksXY(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack); // Compute track matching
   //// Position & Angles
   GlobalMuonTrack matchMFT_MCH_TracksXYPhiTanl(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack);
   //// Position, Angles & Charged Momentum
   GlobalMuonTrack matchMFT_MCH_TracksFull(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack);
   void setMatchingFunction(GlobalMuonTrack (MUONMatching::*func)(GlobalMuonTrack&, MFTTrack&)) { mMatchFunc = func; }
   void setCustomMatchingFunction(GlobalMuonTrack (*func)(GlobalMuonTrack&, MFTTrack&)) { mCustomMatchFunc = func; }
   void runHeavyMatching(); //Finds best match (no search window)
   void fitTracks(); //Fit all matched tracks


 private:

   // Private IO methods
   void loadMFTClusters();

   // Track methods
   bool propagateMCHTrackToMFT(MCHTrack& track); //Propagates MCH Track to Last MFT Plane;
   GlobalMuonTrack MCHtoGlobal(MCHTrack&); // Convert MCH Track to GlobalMuonTrack;



   // Global Muon Track Methods
   void fitGlobalMuonTrack(GlobalMuonTrack&); // Kalman filter
   bool computeCluster(GlobalMuonTrack&, MFTCluster&);
   GlobalMuonTrack (MUONMatching::*mMatchFunc)(GlobalMuonTrack&,MFTTrack&) ;
   GlobalMuonTrack (*mCustomMatchFunc)(GlobalMuonTrack&,MFTTrack&) = nullptr;

   // Data Members
   std::vector<MFTTrack> mMFTTracks;
   std::vector<MCHTrack> mMCHTracks;
   std::vector<MFTTrack> mMCHTracksDummy; // Dummy MCH tracks at the MFT coordinate system
   std::vector<GlobalMuonTrack> mGlobalMuonTracks;
   std::vector<MFTCluster> mMFTClusters;
   std::vector<int> mtrackExtClsIDs;


   o2::dataformats::MCTruthContainer<o2::MCCompLabel> mftTrackLabels;
   o2::dataformats::MCTruthContainer<o2::MCCompLabel> mchTrackLabels;
   o2::dataformats::MCTruthContainer<o2::MCCompLabel> mGlobalTrackLabels;


   // MCH Track Propagation clasee
   o2::mch::TrackExtrap mMCHTrackExtrap;

   double mField_z;
   const double sLastMFTPlaneZ = -77.5;
   double mMatchingPlaneZ = sLastMFTPlaneZ;

};

#include "MUONMatching.cxx"


#endif /* MUON_MATCHING */
