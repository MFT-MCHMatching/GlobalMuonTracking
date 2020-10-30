
#ifndef ALICEO2_GLOBALMUONTRACK_H
#define ALICEO2_GLOBALMUONTRACK_H

#include "Math/SMatrix.h"
#include <TMath.h>
#include <vector>

#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::track {

class GlobalMuonTrack : public o2::track::TrackParCovFwd {
  using ClusRefs = o2::dataformats::RangeRefComp<4>;
  using SMatrix55 =
      ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

public:
  GlobalMuonTrack() = default;
  GlobalMuonTrack(const GlobalMuonTrack &t) = default;
  ~GlobalMuonTrack() = default;

  std::uint32_t getROFrame() const { return mROFrame; }
  void setROFrame(std::uint32_t f) { mROFrame = f; }

  void setMatchingChi2(double chi2) { mMatchingChi2 = chi2; }
  double getMatchingChi2() { return mMatchingChi2; }
  void countCandidate() { mNMFTCandidates++; }
  int getNMFTCandidates() { return mNMFTCandidates; }

  void setBestMFTTrackMatchID(int ID) { mBestMFTTrackMatchID = ID; }
  double getBestMFTTrackMatchID() { return mBestMFTTrackMatchID; }
  void setCloseMatch() { mCloseMatch = true; }
  bool closeMatch() { return mCloseMatch; }

  void print() const;

private:
  std::uint32_t mROFrame = 0; ///< RO Frame
  double mMatchingChi2 = 1.0E308;
  int mBestMFTTrackMatchID = -1;
  int mNMFTCandidates = 0; // Number of candidates within search cut
  bool mCloseMatch = false;
};
} // namespace o2::track

//_________________________________________________________________________________________________
struct MatchingHelper {
  std::string Generator;
  std::string GeneratorConfig;
  std::string MatchingFunction;
  std::string MatchingCutFunc;
  std::string MatchingCutConfig;
  int nMCHTracks = -1;
  int nCorrectMatches = -1;
  int nFakes = -1;
  int nNoMatch = -1;
  int nCloseMatches = -1;
  double matchingPlaneZ;

  int nGMTracks() { return nMCHTracks - nNoMatch; }
  double getCorrectMatchRatio() { return 1.f * nCorrectMatches / nGMTracks(); }
  double getPairingEfficiency() { return 1.f * nGMTracks() / nMCHTracks; }
  std::string Annotation() {
    return Generator + GeneratorConfig + MatchingFunction + "_Z" +
           std::to_string(matchingPlaneZ) + MatchingCutFunc + MatchingCutConfig;
  }
  std::string MatchingConfig() {
    return MatchingFunction + "_Z" + std::to_string(matchingPlaneZ) +
           MatchingCutFunc + MatchingCutConfig;
  }
};

#endif // ALICEO2_GLOBALMUONTRACK_H
