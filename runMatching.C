#if !defined(__CLING__) || defined(__ROOTCLING__)

#endif

#ifdef __MAKECINT__
#pragma link C++ class GlobalMuonTrack+;
#pragma link C++ class std::vector<GlobalMuonTrack>+;
#pragma link C++ class tempMCHTrack+;
#pragma link C++ class std::vector<tempMCHTrack>+;
#endif

#include "MUONMatching.h"

MUONMatching matcher;

//_________________________________________________________________________________________________
// Sample custom matching function that can be passed to MUONMatching
GlobalMuonTrack MyMatchingFunc (GlobalMuonTrack& mchTrack, MFTTrack& mftTrack) {
    auto dx = mchTrack.getX() - mftTrack.getX();
    auto dy = mchTrack.getY() - mftTrack.getY();
    auto score = dx*dx + dy*dy;
    GlobalMuonTrack matchTrack;
    matchTrack.setMatchingChi2(score);
    return matchTrack;
 };


//_________________________________________________________________________________________________
int runMatching()  {


//Custom matching function
//matcher.setCustomMatchingFunction(&MyMatchingFunc);

// Built-in matching functions
//matcher.setMatchingFunction(&MUONMatching::matchMFT_MCH_TracksXY);
//matcher.setMatchingFunction(&MUONMatching::matchMFT_MCH_TracksXYPhiTanl);
matcher.setMatchingFunction(&MUONMatching::matchMFT_MCH_TracksFull);
//matcher.SetMatchingPlane(-45.3);

matcher.loadMFTTracksOut();

//matcher.loadMCHTracks();
//matcher.initGlobalTracks();

matcher.loadDummyMCHTracks();
matcher.initDummyGlobalTracks();

matcher.runHeavyMatching();

matcher.fitTracks();
matcher.saveGlobalMuonTracks();


std::cout << " *** Matching Summary ***" << std::endl;
auto globalTrackID=0;
for (auto gTrack: matcher.getGlobalMuonTracks() ) {
  if (globalTrackID < 15) std::cout << "Best match to MCH Track " << globalTrackID << " is MFT track " << gTrack.getBestMFTTrackMatchID() << " with chi^2 = " <<  gTrack.getMatchingChi2() << std::endl;
  globalTrackID++;
}

return 0;

}
