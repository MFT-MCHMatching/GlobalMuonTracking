
#ifndef TEMP_MCH_ESDTRACKGETTER
#define TEMP_MCH_ESDTRACKGETTER
#include "AliESDMuonTrack.h"
#include "./tempMCHTrack.h"

struct tempMCHTrackGetter : public  AliESDMuonTrack {

  tempMCHTrackGetter(AliESDMuonTrack& t)
        : AliESDMuonTrack{ t } {

         }
  void update(tempMCHTrack& local) {
    local.fInverseBendingMomentum = this->GetInverseBendingMomentum();
    local.fThetaX = this->GetThetaX();
    local.fThetaY = this->GetThetaY();
    local.fZ = this->GetZ();
    local.fBendingCoor = this->GetBendingCoor();
    local.fNonBendingCoor = this->GetNonBendingCoor();
    local.fLabel = this->GetLabel();
    for (int i = 0 ; i < 15 ; i++ ) local.fCovariances[i] = this->fCovariances[i];
  }
  //Double32_t getCovariances() {return  fCovariances;}

};



#endif
