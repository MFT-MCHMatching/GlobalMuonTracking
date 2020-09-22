
#ifndef MUON_MATCHING
#define MUON_MATCHING

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TGeoGlobalMagField.h>
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "MFTTracking/TrackCA.h"
#include "MFTTracking/Cluster.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "TFile.h"
#include "TTree.h"
#include "Math/SMatrix.h"
#endif

#include "include/TrackExtrap.h"
#include "include/GlobalMuonTrack.h"
#include "MCHTracking/TrackParam.h"

using MCHTrack = o2::mch::TrackParam;
using MFTTrack = o2::mft::TrackMFT;
using GlobalMuonTrack = o2::track::GlobalMuonTrack;
using MCLabels = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
using MFTClusters = o2::mft::Cluster;

class MUONMatching
{
 public:
   MUONMatching();
   ~MUONMatching() = default;
   void Clear();
   void SetMatchingPlane(double z) {sMatchingPlaneZ = z;}

   // Track IO
   void loadMCHTracks();
   void loadDummyMCHTracks();
   void loadMFTTracksOut();
   void saveGlobalMuonTracks();
   std::vector<GlobalMuonTrack> getGlobalMuonTracks() const { return mGlobalMuonTracks;}

   void loadROFrameData(int); //Loads data from
   void initGlobalTracks(); //


   // Matching methods
   void runHeavyMatching(); //Finds best match (no search window)
   void fitTracks(); //Fit all matched tracks


 private:

   // Track methods
   bool propagateMCHTrackToMFT(MCHTrack& track); //Propagates MCH Track to Last MFT Plane;
   GlobalMuonTrack MCHtoGlobal(MCHTrack&); // Convert MCH Track to GlobalMuonTrack;

   // Matching methods
   // Position
   GlobalMuonTrack matchMFT_MCH_TracksXY(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack); // Compute track matching; returns matching chi2
   double matchMFT_MCH_TracksXY(MCHTrack& mchTrack, MFTTrack& mftTrack); // Compute track matching; returns matching chi2
   //// Position & Angles
   double matchMFT_MCH_TracksXYPhiTanl(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack); // Compute track matching; returns matching chi2
   double matchMFT_MCH_TracksXYPhiTanl(MCHTrack& mchTrack, MFTTrack& mftTrack); // Compute track matching; returns matching chi2
   //// Position, Angles & Charged Momentum
   double matchMFT_MCH_TracksFull(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack); // Compute track matching; returns matching chi2
   double matchMFT_MCH_TracksFull(MCHTrack& mchTrack, MFTTrack& mftTrack); // Compute track matching; returns matching chi2

   // Global Muon Track Methods
   void fitGlobalMuonTrack(GlobalMuonTrack&); // Kalman filter


   // Data Members
   std::vector<MFTTrack> mMFTTracks;
   std::vector<MCHTrack> mMCHTracks;
   std::vector<MFTTrack> mMCHTracksDummy; // Dummy MCH at the MFT coordinate system
   std::vector<GlobalMuonTrack> mGlobalMuonTracks;
   std::vector<MFTClusters> mMFTClusters;


   o2::dataformats::MCTruthContainer<o2::MCCompLabel> mftTrackLabels;
   o2::dataformats::MCTruthContainer<o2::MCCompLabel> mchTrackLabels;

   // MCH Track Propagation clasee
   o2::mch::TrackExtrap mMCHTrackExtrap;

   double mField_z;
   const double sLastMFTPlaneZ = -77.5;
   double sMatchingPlaneZ = sLastMFTPlaneZ;

};

#include "MUONMatching.cxx"


#endif /* MUON_MATCHING */
