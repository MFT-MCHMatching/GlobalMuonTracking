#ifndef AliGenMimicPbPb_H
#define AliGenMimicPbPb_H

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

#include <fstream>
#include <iostream>

class AliGenMimicPbPb : public AliGenerator {

public:
  
  AliGenMimicPbPb();
  AliGenMimicPbPb(Int_t nPart/*, Char_t *inputFile*/);

  virtual ~AliGenMimicPbPb() {}
  virtual void Generate();
  virtual void Init();
  
private:

  AliGenMimicPbPb(const AliGenMimicPbPb&);
  AliGenMimicPbPb &operator=(const AliGenMimicPbPb&);

protected:

  TH1F* fHistNumFwdPrimePion;
  TH1F* fHistNumFwdPrimeKaon;
  TH1F* fHistNumFwdPrimeProton;
  TH1F* fHistNumFwdPrimeMuon;
  TH1F* fHistNumFwdPrimeElectron;

  TH2F* fHistPrimePionRapPt;
  TH2F* fHistPrimeKaonRapPt;
  TH2F* fHistPrimeProtonRapPt;
  TH2F* fHistPrimeMuonRapPt;
  TH2F* fHistPrimeElectronRapPt;
  
  Int_t fNType;

  TH1F* fNum[10];
  TH2F* fDist[10];

  Int_t fPdgCode[10];
  Int_t fNPart[10];
  Double_t fMass[10];
  
  ClassDef(AliGenMimicPbPb, 1)

};

//====================================================================================================================================================

#endif


