
#ifndef ALICEO2_GLOBALMUONTRACK_H
#define ALICEO2_GLOBALMUONTRACK_H

#include "Math/SMatrix.h"
#include <TMath.h>
#include <vector>

#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "SimulationDataFormat/MCCompLabel.h"

//_________________________________________________________________________________________________
struct MatchingHelper {
  std::string Generator;
  std::string GeneratorConfig;
  std::string MatchingFunction;
  std::string MatchingCutFunc;
  std::string MatchingCutConfig;
  std::string MLFeaturesFunction;
  int nMCHTracks = -1;
  int nCorrectMatches = -1;
  int nFakes = -1;
  int nNoMatch = -1;
  int nCloseMatches = -1;
  double matchingPlaneZ;

  int nGMTracks() { return nMCHTracks - nNoMatch; }
  double getCorrectMatchRatio() { return 1.f * nCorrectMatches / nGMTracks(); }
  double getPairingEfficiency() { return 1.f * nGMTracks() / nMCHTracks; }
  std::string Annotation()
  {
    return Generator + GeneratorConfig + MatchingFunction + (MLFeaturesFunction == "" ? "" : "_") + MLFeaturesFunction + "_Z" +
           std::to_string(matchingPlaneZ) + MatchingCutFunc + MatchingCutConfig;
  }
  std::string MatchingConfig()
  {
    return MatchingFunction + (MLFeaturesFunction == "" ? "" : "_") + MLFeaturesFunction + "_Z" + std::to_string(matchingPlaneZ) +
           MatchingCutFunc + MatchingCutConfig;
  }
};

#endif // ALICEO2_GLOBALMUONTRACK_H
