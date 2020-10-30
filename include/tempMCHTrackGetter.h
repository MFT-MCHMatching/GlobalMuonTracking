
#ifndef TEMP_MCH_ESDTRACKGETTER
#define TEMP_MCH_ESDTRACKGETTER
#include "./tempMCHTrack.h"
#include "AliESDMuonTrack.h"

struct tempMCHTrackGetter : public AliESDMuonTrack {

  tempMCHTrackGetter(AliESDMuonTrack &t) : AliESDMuonTrack{t} {}
  void update(tempMCHTrack &local) {
    local.fInverseBendingMomentum =
        this->GetInverseBendingMomentumUncorrected();
    local.fThetaX = this->GetThetaXUncorrected();
    local.fThetaY = this->GetThetaYUncorrected();
    local.fZ = this->GetZUncorrected();
    local.fBendingCoor = this->GetBendingCoorUncorrected();
    local.fNonBendingCoor = this->GetNonBendingCoorUncorrected();
    local.fLabel = this->GetLabel();
    for (int i = 0; i < 15; i++)
      local.fCovariances[i] = this->fCovariances[i];
  }
  // Double32_t getCovariances() {return  fCovariances;}
};

#endif
