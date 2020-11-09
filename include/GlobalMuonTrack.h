
#ifndef ALICEO2_GLOBALMUONTRACK_H
#define ALICEO2_GLOBALMUONTRACK_H

#include "Math/SMatrix.h"
#include <TMath.h>
#include <vector>

#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::track
{

class GlobalMuonTrack : public o2::track::TrackParCovFwd
{
  using ClusRefs = o2::dataformats::RangeRefComp<4>;
  using SMatrix55 =
    ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

 public:
  GlobalMuonTrack() = default;
  GlobalMuonTrack(const GlobalMuonTrack& t) = default;
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
  void computeResiduals2Cov(const o2::track::TrackParCovFwd& t)
  {
    mResiduals2Cov(0) =
      (getX() - t.getX()) /
      TMath::Sqrt(getCovariances()(0, 0) + t.getCovariances()(0, 0));
    mResiduals2Cov(1) =
      (getY() - t.getY()) /
      TMath::Sqrt(getCovariances()(1, 1) + t.getCovariances()(1, 1));
    mResiduals2Cov(2) =
      (getPhi() - t.getPhi()) /
      TMath::Sqrt(getCovariances()(2, 2) + t.getCovariances()(2, 2));
    mResiduals2Cov(3) =
      (getTanl() - t.getTanl()) /
      TMath::Sqrt(getCovariances()(3, 3) + t.getCovariances()(3, 3));
    mResiduals2Cov(4) =
      (getInvQPt() - t.getInvQPt()) /
      TMath::Sqrt(getCovariances()(4, 4) + t.getCovariances()(4, 4));
    ;
  }
  const SMatrix5& getResiduals2Cov() { return mResiduals2Cov; }

  void print() const;

  /// for performance studies: store parameters and covariances for the MCH and MFT
  const SMatrix5& getParametersMCH() const { return mParametersMCH; }
  const SMatrix5& getParametersMFT() const { return mParametersMFT; }
  void setParametersMCH(const SMatrix5& parameters) { mParametersMCH = parameters; }
  void setParametersMFT(const SMatrix5& parameters) { mParametersMFT = parameters; }
  const SMatrix55& getCovariancesMCH() const { return mCovariancesMCH; }
  const SMatrix55& getCovariancesMFT() const { return mCovariancesMFT; }
  void setCovariancesMCH(const SMatrix55& covariances)
  {
    mCovariancesMCH = covariances;
  }
  void setCovariancesMFT(const SMatrix55& covariances)
  {
    mCovariancesMFT = covariances;
  }
  void setMCHTrackID(int ID) { mMCHTrackID = ID; }
  double getMCHTrackID() { return mMCHTrackID; }

 private:
  std::uint32_t mROFrame = 0; ///< RO Frame
  double mMatchingChi2 = 1.0E308;
  int mBestMFTTrackMatchID = -1;
  int mNMFTCandidates = 0; // Number of candidates within search cut
  bool mCloseMatch = false;

  SMatrix5 mResiduals2Cov;

  /// for performance studies: store parameters and covariances for the MCH and MFT
  SMatrix5 mParametersMCH{};   ///< \brief Track parameters, MCH part of the track
  SMatrix5 mParametersMFT{};   ///< \brief Track parameters, MFT part of the track
  SMatrix55 mCovariancesMCH{}; ///< \brief Covariance matrix of track parameters, MCH
  SMatrix55 mCovariancesMFT{}; ///< \brief Covariance matrix of track parameters, MFT
  int mMCHTrackID = -1;
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
  std::string Annotation()
  {
    return Generator + GeneratorConfig + MatchingFunction + "_Z" +
           std::to_string(matchingPlaneZ) + MatchingCutFunc + MatchingCutConfig;
  }
  std::string MatchingConfig()
  {
    return MatchingFunction + "_Z" + std::to_string(matchingPlaneZ) +
           MatchingCutFunc + MatchingCutConfig;
  }
};

#endif // ALICEO2_GLOBALMUONTRACK_H
