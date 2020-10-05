#if !defined(__CLING__) || defined(__ROOTCLING__)

#endif

#include "MUONMatcher.h"

//#ifdef __MAKECINT__
#pragma link C++ class GlobalMuonTrack+;
#pragma link C++ class std::vector<GlobalMuonTrack>+;
#pragma link C++ class tempMCHTrack+;
#pragma link C++ class std::vector<tempMCHTrack>+;
//#endif


MUONMatcher matcher;

//_________________________________________________________________________________________________
// Sample custom matching function that can be passed to MUONMatcher
double MyMatchingFunc (const GlobalMuonTrack& mchTrack, const MFTTrack& mftTrack) {
    auto dx = mchTrack.getX() - mftTrack.getX();
    auto dy = mchTrack.getY() - mftTrack.getY();
    auto score = dx*dx + dy*dy;
    return score;
 };


 //_________________________________________________________________________________________________
 // Sample custom cut criteria that can be passed to MUONMatcher
bool MyMatchingCut (GlobalMuonTrack& mchTrack, MFTTrack& mftTrack) {
   auto cutDistance = 1.0;
   auto dx = mchTrack.getX() - mftTrack.getX();
   auto dy = mchTrack.getY() - mftTrack.getY();
   auto distance = TMath::Sqrt(dx*dx + dy*dy);
   return distance < cutDistance;
  };

//_________________________________________________________________________________________________
int runMatching()  {

gSystem->Load("libO2MCHTracking");
//Custom matching function
//matcher.setCustomMatchingFunction(&MyMatchingFunc, "_aliasForMyMatchingFunction");

// Built-in matching functions
//matcher.setMatchingFunction(&MUONMatcher::matchMFT_MCH_TracksXY);
//matcher.setMatchingFunction(&MUONMatcher::matchMFT_MCH_TracksXYPhiTanl);
matcher.setMatchingFunction(&MUONMatcher::matchMFT_MCH_TracksFull);
//matcher.SetMatchingPlaneZ(-45.3);
//matcher.SetMatchingPlaneZ(0.);
//
//matcher.setCutFunction(&MUONMatcher::matchCutDistance);
//matcher.setCutDistanceParam(1.0);
//
//
//matcher.setCutFunction(&MUONMatcher::matchCutDistanceSigma);
//matcher.setCutDistanceParam(3.0);

matcher.SetMatchingPlaneZ(-80.0);
//matcher.SetVerbosity(true);


matcher.loadMFTTracksOut(); // Load and propagates to matching plane
matcher.loadMCHTracks();
matcher.initGlobalTracks(); // Propagate MCH tracks to matching plane and convert parameters and covariances matrix to MFT coordinate system
//matcher.loadDummyMCHTracks();
//matcher.initDummyGlobalTracks();
matcher.runEventMatching(); // Runs track matching event-by-event
matcher.fitTracks(); // Kalman filter
matcher.saveGlobalMuonTracks();

/*
std::cout << " *** Matching Summary ***" << std::endl;
auto globalTrackID=0;
for (auto gTrack: matcher.getGlobalMuonTracks() ) {
  if (globalTrackID < 15) std::cout << "Best match to MCH Track " << globalTrackID << " is MFT track " << gTrack.getBestMFTTrackMatchID() << " with chi^2 = " <<  gTrack.getMatchingChi2() << std::endl;
  globalTrackID++;
}
*/
return 0;

}
