#if !defined(__CLING__) || defined(__ROOTCLING__)


#endif

#ifdef __MAKECINT__
#pragma link C++ class GlobalMuonTrack+;
#pragma link C++ class std::vector<GlobalMuonTrack>+;
#endif


#include "MUONMatching.h"

MUONMatching matcher;

int runMatching()  {

matcher.loadMFTTracksOut();
//matcher.loadMCHTracks();
matcher.loadDummyMCHTracks();
matcher.initGlobalTracks();
matcher.runHeavyMatching();
matcher.saveGlobalMuonTracks();

std::cout << " *** Matching Summary ***" << std::endl;
auto globalTrackID=0;
for (auto gTrack: matcher.getGlobalMuonTracks() ) {
  if (globalTrackID < 15) std::cout << "Best match to MCH Track " << globalTrackID << " is MFT track " << gTrack.getBestMFTTrackMatchID() << " with chi^2 = " <<  gTrack.getMatchingChi2() << std::endl;
  globalTrackID++;
}

return 0;

}
