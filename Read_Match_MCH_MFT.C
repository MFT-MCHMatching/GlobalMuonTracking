#if !defined(__CLING__) || defined(__ROOTCLING__)

#ifdef __MAKECINT__
#pragma link C++ class GlobalMuonTrack + ;
#pragma link C++ class std::vector < GlobalMuonTrack> + ;
#pragma link C++ class GlobalMuonTrackExt + ;
#pragma link C++ class std::vector < GlobalMuonTrackExt> + ;
#pragma link C++ class MatchingHelper + ;
#endif

#include <TFile.h>
#include <TTree.h>
#include <Math/SMatrix.h>
#include "include/GlobalMuonTrack.h"

#endif

using GlobalMuonTrack = o2::track::GlobalMuonTrack;
using GlobalMuonTrackExt = o2::track::GlobalMuonTrackExt;
using SMatrix5 = o2::track::SMatrix5;
using SMatrix55Sym = o2::track::SMatrix55Sym;

void Read_Match_MCH_MFT(
  const std::string trkFile = "GlobalMuonTracksExt.root")
{

  // Global Muon Tracks
  TFile* trkFileIn = new TFile(trkFile.c_str());
  TTree* gmTrackTree = (TTree*)trkFileIn->Get("o2sim");
  std::vector<GlobalMuonTrackExt> trackGMVec, *trackGMVecP = &trackGMVec;
  gmTrackTree->SetBranchAddress("GlobalMuonTrackExt", &trackGMVecP);

  Int_t numberOfEvents = gmTrackTree->GetEntries();
  Int_t igm = 0;
  for (int iEvent = 0; iEvent < numberOfEvents; iEvent++) {

    gmTrackTree->GetEntry(iEvent);
    printf("\n\n        ##### Event %i #####\n\n", iEvent);
    for (auto& gmTrack : trackGMVec) {
      const SMatrix5& trackMCHpar = gmTrack.getParametersMCH();
      const SMatrix5& trackMFTpar = gmTrack.getParametersMFT();
      const SMatrix55Sym& trackMCHcov = gmTrack.getCovariancesMCH();
      const SMatrix55Sym& trackMFTcov = gmTrack.getCovariancesMFT();
      if (gmTrack.closeMatch()) {
        printf("Global track %d close match!\n", igm);
      } else {
        printf("Global track %d \n", igm);
      }
      printf("MCH parameters (x, y, phi, tanl, invqpt):\n");
      for (auto& par : trackMCHpar) {
        printf("%f ", par);
      }
      printf("\n");
      printf("MFT parameters (x, y, phi, tanl, invqpt):\n");
      for (auto& par : trackMFTpar) {
        printf("%f ", par);
      }
      printf("\n");
      printf("MCH covariances:\n");
      for (Int_t i = 0; i < 5; i++) {
        for (Int_t j = 0; j < 5; j++) {
          printf("%f ", trackMCHcov(i, j));
        }
        printf("\n");
      }
      printf("MFT covariances:\n");
      for (Int_t i = 0; i < 5; i++) {
        for (Int_t j = 0; j < 5; j++) {
          printf("%f ", trackMFTcov(i, j));
        }
        printf("\n");
      }
      printf("-------------------------------------------------------------\n");
      ++igm;
    }
  }
  trkFileIn->Close();
}
