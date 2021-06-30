#ifndef AliGenDimuon_H
#define AliGenDimuon_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Parametric generator of uncorrelated muons
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliGenerator.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

class AliGenDimuon : public AliGenerator {

public:
  
  AliGenDimuon();
  AliGenDimuon(Int_t nPart/*, Char_t *inputFile*/);

  virtual ~AliGenDimuon() {}
  virtual void Generate();
  virtual void Init();

  void SetGenerateParticle(Int_t code){fPdgCode = code;}
  
  Bool_t SetMuonMomentumRange(Double_t min, Double_t max){
    kMuonMomentumRange = kTRUE;
    fMinMuonMomentum   = min;
    fMaxMuonMomentum   = max;
  }

  Bool_t SetMuonRapRange(Double_t min, Double_t max){
    kMuonRapRange = kTRUE;
    fMinMuonRap   = min;
    fMaxMuonRap   = max;
  }

  Bool_t SetMuonEtaRange(Double_t min, Double_t max){
    kMuonEtaRange = kTRUE;
    fMinMuonEta   = min;
    fMaxMuonEta   = max;
  }

private:

  AliGenDimuon(const AliGenDimuon&);
  AliGenDimuon &operator=(const AliGenDimuon&);

  Int_t fPdgCode;
  
  Double_t fMinMuonMomentum;
  Double_t fMaxMuonMomentum;

  Double_t fMinMuonRap;
  Double_t fMaxMuonRap;

  Double_t fMinMuonEta;
  Double_t fMaxMuonEta;
  
  Bool_t kMuonMomentumRange;
  Bool_t kMuonRapRange;
  Bool_t kMuonEtaRange;
  
  TH2F* fHistParentKine;
  
  Double_t EtaToTheta(Double_t arg);

protected:

  //enum E_Parent{Rho2Body,Omega2Body,Phi2Body,Jpsi2Body,Psiprime,Upsilon1S,Upsilon2S,Upsilon3S,BtoJpasi,LMRinMedium};

  Double_t mass[10];
  
  ClassDef(AliGenDimuon, 1)

};

//====================================================================================================================================================

#endif


