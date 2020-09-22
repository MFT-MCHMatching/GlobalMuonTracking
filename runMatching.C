#if !defined(__CLING__) || defined(__ROOTCLING__)


#endif

#include "MUONMatching.h"

MUONMatching matcher;

int runMatching()  {

matcher.loadMFTTracksOut();
//matcher.loadMCHTracks();
matcher.loadDummyMCHTracks();
matcher.initGlobalTracks();
matcher.runHeavyMatching();
return 0;

}
