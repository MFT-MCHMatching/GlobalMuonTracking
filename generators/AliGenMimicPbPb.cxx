/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//====================================================================================================================================================
//
//      
//
//      
//
//====================================================================================================================================================

#include "AliConst.h"
#include "AliRun.h"
#include "AliGenEventHeader.h"
#include "TDatabasePDG.h"
#include "AliPDG.h"
#include "TFile.h"
#include "TF1.h"
#include "TROOT.h"
#include "AliGenMimicPbPb.h"
#include "TVector3.h"
#include "AliLog.h"
#include "TGenPhaseSpace.h"

#include <fstream>
#include <iostream>
#include <unistd.h>

using namespace std ;

ClassImp(AliGenMimicPbPb)

//====================================================================================================================================================

AliGenMimicPbPb::AliGenMimicPbPb():
  AliGenerator(), 
  fHistNumFwdPrimePion(NULL),
  fHistNumFwdPrimeKaon(NULL),
  fHistNumFwdPrimeProton(NULL),
  fHistNumFwdPrimeMuon(NULL),
  fHistNumFwdPrimeElectron(NULL),
  fHistPrimePionRapPt(NULL),
  fHistPrimeKaonRapPt(NULL),
  fHistPrimeProtonRapPt(NULL),
  fHistPrimeMuonRapPt(NULL),
  fHistPrimeElectronRapPt(NULL),
  fNType(5),
  fDist(),
  fNum(),
  fPdgCode(),
  fNPart(),
  fMass()
{
  // Default constructor    
}

//====================================================================================================================================================

AliGenMimicPbPb::AliGenMimicPbPb(Int_t nPart/*, Char_t *inputFile*/):
  AliGenerator(nPart),
  fHistNumFwdPrimePion(NULL),
  fHistNumFwdPrimeKaon(NULL),
  fHistNumFwdPrimeProton(NULL),
  fHistNumFwdPrimeMuon(NULL),
  fHistNumFwdPrimeElectron(NULL),
  fHistPrimePionRapPt(NULL),
  fHistPrimeKaonRapPt(NULL),
  fHistPrimeProtonRapPt(NULL),
  fHistPrimeMuonRapPt(NULL),
  fHistPrimeElectronRapPt(NULL),
  fNType(5),
  fDist(),
  fNum(),
  fPdgCode(),
  fNPart(),
  fMass()
{

  Init();
  // Standard constructor
 
}

//====================================================================================================================================================

void AliGenMimicPbPb::Generate() {

  // Generate one trigger
  
  Double_t polar[3]= {0,0,0};
  Int_t nt;
  Double_t origin[3];
  Double_t time=0;
  
  for (Int_t j=0; j<3; j++) origin[j] = fOrigin[j];
  time = fTimeOrigin;
  if (fVertexSmear==kPerEvent) {
    Vertex();
    for (Int_t j=0; j<3; j++) origin[j] = fVertex[j];
    time = fTime;
  }
  
  Int_t nPartGenerated = 0;
  
  TLorentzVector particle;
  
  for(Int_t iType=0; iType<fNType; ++iType){
    
    fNPart[iType] = fNum[iType] -> GetRandom();

    for(Int_t iPart=0; iPart<fNPart[iType]; ++iPart){
      
      Double_t pt   = 0;
      Double_t eta  = 0;
      Double_t phi  = gRandom->Uniform(0.,TMath::TwoPi());
      
      fDist[iType]->GetRandom2(eta,pt);
      
      Int_t charge = 1;
      
      if (gRandom->Rndm() < 0.5){
	charge = +1;
      }
      else{
	charge = -1;
      }

      particle.SetPtEtaPhiM(pt,eta,phi,fMass[iType]);
      
      Double_t theta = particle.Theta();
      
      if (TestBit(kThetaRange) && (theta<fThetaMin || theta>fThetaMax)){
	continue;
      }
      
      PushTrack(1, -1, charge * fPdgCode[iType],
		particle.Px(),particle.Py(),particle.Pz(),particle.E(),
		origin[0],origin[1],origin[2],Double_t(time),
		polar[0],polar[1],polar[2],
		kPPrimary, nt, 1., 1);
      
      nPartGenerated++;
      
    }

  }
  
  AliGenEventHeader* header = new AliGenEventHeader("Mimic_HIJING");
  header->SetPrimaryVertex(fVertex);
  header->SetNProduced(nPartGenerated);
  header->SetInteractionTime(fTime);
  
  // Passes header either to the container or to gAlice
  if (fContainer) {
    fContainer->AddHeader(header);
  } 
  else {
    gAlice->SetGenEventHeader(header);	
  }
  
}

//====================================================================================================================================================

void AliGenMimicPbPb::Init() {

  // Initialisation, check consistency of selected ranges
  /*
  if (TestBit(kPtRange) && TestBit(kMomentumRange)) 
    Fatal("Init","You should not set the momentum range and the pt range at the same time!\n");
  if ((!TestBit(kPtRange)) && (!TestBit(kMomentumRange))) 
    Fatal("Init","You should set either the momentum or the pt range!\n");
  if ((TestBit(kYRange) && TestBit(kThetaRange)) || (TestBit(kYRange) && TestBit(kEtaRange)) || (TestBit(kEtaRange) && TestBit(kThetaRange)))
    Fatal("Init","You should only set the range of one of these variables: y, eta or theta\n");
  if ((!TestBit(kYRange)) && (!TestBit(kEtaRange)) && (!TestBit(kThetaRange)))
    Fatal("Init","You should set the range of one of these variables: y, eta or theta\n");
  */

  AliPDG::AddParticlesToPdgDataBase();

  TFile* inFile = TFile::Open("./inputHijingParam.root");

  fNum[0] = (TH1F*)inFile->Get("fHistNumFwdPrimePion")     -> Clone();
  fNum[1] = (TH1F*)inFile->Get("fHistNumFwdPrimeKaon")     -> Clone();
  fNum[2] = (TH1F*)inFile->Get("fHistNumFwdPrimeProton")   -> Clone();
  fNum[3] = (TH1F*)inFile->Get("fHistNumFwdPrimeMuon")     -> Clone();
  fNum[4] = (TH1F*)inFile->Get("fHistNumFwdPrimeElectron") -> Clone();
    
  fDist[0] = (TH2F*)inFile->Get("fHistPrimePionRapPt")     -> Clone();
  fDist[1] = (TH2F*)inFile->Get("fHistPrimeKaonRapPt")     -> Clone();
  fDist[2] = (TH2F*)inFile->Get("fHistPrimeProtonRapPt")   -> Clone();
  fDist[3] = (TH2F*)inFile->Get("fHistPrimeMuonRapPt")     -> Clone();
  fDist[4] = (TH2F*)inFile->Get("fHistPrimeElectronRapPt") -> Clone();
  
  fMass[0] = TDatabasePDG::Instance()->GetParticle(211)  -> Mass();
  fMass[1] = TDatabasePDG::Instance()->GetParticle(321)  -> Mass();
  fMass[2] = TDatabasePDG::Instance()->GetParticle(2212) -> Mass();
  fMass[3] = TDatabasePDG::Instance()->GetParticle(13)   -> Mass();
  fMass[4] = TDatabasePDG::Instance()->GetParticle(11)   -> Mass();

  fPdgCode[0] = 211;
  fPdgCode[1] = 321;
  fPdgCode[2] = 2212;
  fPdgCode[3] = 13;
  fPdgCode[4] = 11;
  
}
