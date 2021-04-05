
#ifndef TEMP_MCH_ESDTRACK
#define TEMP_MCH_ESDTRACK

struct tempMCHTrack {
  // A transitional class to get MCH tracks from aliroot to O2

  Double32_t fInverseBendingMomentum;
  Double32_t fThetaX;
  Double32_t fThetaY;
  Double32_t fZ;
  Double32_t fBendingCoor;
  Double32_t fNonBendingCoor;
  Double32_t fCovariances[15];
  Int_t fLabel;
  Int_t fiEv;
};

#endif
