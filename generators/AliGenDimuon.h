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
  
private:

  AliGenDimuon(const AliGenDimuon&);
  AliGenDimuon &operator=(const AliGenDimuon&);

  void SetPtShape();
  void SetRapidityShape();

  Int_t fPdgCode;

protected:

  TF1* fPt;
  TF1* fRap;

  //enum E_Parent{Rho2Body,Omega2Body,Phi2Body,Jpsi2Body,Psiprime,Upsilon1S,Upsilon2S,Upsilon3S,BtoJpasi,LMRinMedium};

  Double_t mass[10];

  TF1* fBW;
  
  ClassDef(AliGenDimuon, 1)

};

//====================================================================================================================================================

#endif


