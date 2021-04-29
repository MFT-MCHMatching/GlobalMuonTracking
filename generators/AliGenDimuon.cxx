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
#include "AliGenDimuon.h"
#include "TVector3.h"
#include "AliLog.h"
#include "TGenPhaseSpace.h"

#include <iostream>

using namespace std ;

ClassImp(AliGenDimuon)

  //====================================================================================================================================================

  AliGenDimuon::AliGenDimuon() : AliGenerator(),
  fPt(0x0),
  fRap(0x0),
  fPdgCode(0),
  fBW(0x0),
  fMinMuonMomentum(4.0),
  fMaxMuonMomentum(9999.),
  fMinMuonRap(-3.8),
  fMaxMuonRap(-2.3),
  fMinMuonEta(-3.8),
  fMaxMuonEta(-2.3),
  kMuonMomentumRange(kTRUE),
  kMuonRapRange(kFALSE),
  kMuonEtaRange(kTRUE),
  fHistParentKine(0x0)
{

  // Default constructor
}

//====================================================================================================================================================

AliGenDimuon::AliGenDimuon(Int_t nPart/*, Char_t *inputFile*/):
  AliGenerator(nPart),
  fPt(0x0),
  fRap(0x0),
  fPdgCode(0),
  fMinMuonMomentum(4.0),
  fMaxMuonMomentum(9999.),
  fMinMuonRap(-3.8),
  fMaxMuonRap(-2.3),
  fMinMuonEta(-3.8),
  fMaxMuonEta(-2.3),
  kMuonMomentumRange(kTRUE),
  kMuonRapRange(kFALSE),
  kMuonEtaRange(kTRUE),
  fHistParentKine(0x0)
{

  // Standard constructor

  fName  = "ParamDimuons";
  fTitle = "Parametric muon pair generator";

  SetPtShape();
  SetRapidityShape();
  
  fBW = new TF1("fBW","gaus",0.2,15);
}

//====================================================================================================================================================

void AliGenDimuon::Generate() {
  
  // Generate one trigger
  Double_t polar[3]= {0,0,0};
  Int_t nt;
  Double_t origin[3];

  Double_t mass=0.;
  Double_t pt=0.;
  Double_t px=0.;
  Double_t py=0.;
  Double_t pz=0.;
  Double_t p=0.;
  Double_t rap=0.;
  Double_t mom=0.;
  Double_t energy=0;
  Double_t phi=0.;
  Double_t time=0.;
  Double_t theta = 0.;
  Double_t eta=0.;

  Int_t pdgCode1;
  Int_t pdgCode2;

  Double_t energy1 = 0;
  Double_t px1     = 0;
  Double_t py1     = 0;
  Double_t pz1     = 0;
  Double_t pt1     = 0;
  Double_t mom1    = 0;
  Double_t rap1    = 0;
  Double_t theta1  = 0;
  Double_t eta1    = 0;

  Double_t energy2 = 0;
  Double_t px2     = 0;
  Double_t py2     = 0;
  Double_t pz2     = 0;
  Double_t pt2     = 0;
  Double_t mom2    = 0;
  Double_t rap2    = 0;
  Double_t theta2  = 0;
  Double_t eta2    = 0;

  for (Int_t j=0; j<3; j++) origin[j] = fOrigin[j];
  time = fTimeOrigin;
  if (fVertexSmear==kPerEvent) {
    Vertex();
    for (Int_t j=0; j<3; j++) origin[j] = fVertex[j];
    time = fTime;
  }

  Int_t nPartGenerated = 2;

  Double_t mass_muon    = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  Double_t mass_parent  = TDatabasePDG::Instance()->GetParticle(fPdgCode)->Mass();
  Double_t width_parent = TDatabasePDG::Instance()->GetParticle(fPdgCode)->Width();

  mass = mass_parent;

  Double_t daughter_m[2]={mass_muon,mass_muon};

  TLorentzVector rest_p(0,0,0,mass);

  TGenPhaseSpace ps_decay;
  ps_decay.SetDecay(rest_p,2,daughter_m);
  
  Int_t nParentTrial = 0;
  Int_t nDaughterTrial = 0;
  
  if (kMuonMomentumRange){
    printf("Set kMuonMomentumRange:  Range(%.2f,%.2f) \n",fMinMuonMomentum,fMaxMuonMomentum);
  }
  if (kMuonRapRange){
    printf("Set kMuonRapRange:  Range(%.2f,%.2f) \n",fMinMuonRap,fMaxMuonRap);
  }
  if (kMuonEtaRange){
    printf("Set kMuonEtaRange:  Range(%.2f,%.2f) \n",fMinMuonEta,fMaxMuonEta);
  }
  

  while (1){

    ++nParentTrial;
    
    /////////////////////////////////////////////////////////////////////////////////
    //Parent particle kinematics
    /////////////////////////////////////////////////////////////////////////////////
    
    //Get parent kinemtacs from histogram
    fHistParentKine->GetRandom2(rap,pt);
    rap = -1*fabs(rap);

    //pt  = fPt  -> GetRandom(0.,1.);
    //rap = fRap -> GetRandom(-3.6,-2.5);
    phi = gRandom->Uniform(0.,TMath::TwoPi());

    px = pt*TMath::Cos(phi);
    py = pt*TMath::Sin(phi);
    energy   = sqrt((px*px + py*py + mass*mass)/(1-pow((exp(2*rap)-1)/(1+exp(2*rap)),2)));
    pz    = (exp(2*rap)-1)/(1+exp(2*rap)) * energy;
    p     = sqrt(px*px + py*py + pz*pz);
    theta = TMath::ACos(pz/p);
    eta = -1*TMath::Log(TMath::Tan(theta/2.));

    TLorentzVector parent;
    parent.SetPtEtaPhiM(pt,eta,phi,mass);

    mom   = parent.P();
    theta = parent.Theta();
    
    //printf("Set kMuonMomentumRange:  Range(%.2f,%.2f) \n",fMinMuonMomentum,fMaxMuonMomentum);
    //printf("Set kMuonRapRange:  Range(%.2f,%.2f) \n",fMinMuonRap,fMaxMuonRap);
    //printf("Set kMuonEtaRange:  Range(%.2f,%.2f) \n",fMinMuonEta,fMaxMuonEta);
    //printf("Generated parent particle kinematics: Pt = %.3f [GeV/c]  y = %.3f\n",pt,rap);

    //Select parent particle kinematics 
    if (TestBit(kPtRange)    && (pt<fPtMin || pt>fPtMax))             continue;
    if (TestBit(kYRange)     && (rap<fYMin || rap>fYMax))             continue;
    
    TVector3 vec_boost =parent.BoostVector();
    
    while (1){

      ++nDaughterTrial;
      
      /////////////////////////////////////////////////////////////////////////////////
      //Daughter muon kinematics
      /////////////////////////////////////////////////////////////////////////////////
      
      ps_decay.Generate();

      TLorentzVector* muon1 = ps_decay.GetDecay(0);
      TLorentzVector* muon2 = ps_decay.GetDecay(1);

      TLorentzVector* boostedMuon1 = (TLorentzVector*)muon1->Clone();
      TLorentzVector* boostedMuon2 = (TLorentzVector*)muon2->Clone();

      //Boost the daughter muons with the parent particle kinematics
      boostedMuon1->Boost(vec_boost);
      boostedMuon2->Boost(vec_boost);

      //Daughter muon1 kinematics in lab frame
      energy1 = boostedMuon1->E();
      px1     = boostedMuon1->Px();
      py1     = boostedMuon1->Py();
      pz1     = boostedMuon1->Pz();
      pt1     = boostedMuon1->Pt();
      mom1    = boostedMuon1->P();
      rap1    = boostedMuon1->Rapidity();
      theta1  = boostedMuon1->Theta();
      eta1    = boostedMuon1->Eta();
      

      //Select daughter muons kinematics
      if (kMuonMomentumRange && (fMinMuonMomentum>mom1 || mom1>fMaxMuonMomentum)){
	delete boostedMuon1;
	delete boostedMuon2;
	continue; 
      }
      if (kMuonRapRange      && (fMaxMuonRap<rap1      || rap1<fMinMuonRap)){
	delete boostedMuon1;
	delete boostedMuon2;
	continue;
      }
      if (kMuonEtaRange      && (fMaxMuonEta<eta1      || eta1<fMinMuonEta)){
	delete boostedMuon1;
	delete boostedMuon2;
	continue;
      }
      
      //Daughter boostedMuon2 kinematics in lab frame
      energy2 = boostedMuon2->E();
      px2     = boostedMuon2->Px();
      py2     = boostedMuon2->Py();
      pz2     = boostedMuon2->Pz();
      pt2     = boostedMuon2->Pt();
      mom2    = boostedMuon2->P();
      rap2    = boostedMuon2->Rapidity();
      theta2  = boostedMuon2->Theta();
      eta2    = boostedMuon1->Eta();

      //Select daughter muons kinematics
      if (kMuonMomentumRange && (fMinMuonMomentum>mom2 || mom2>fMaxMuonMomentum)){
	delete boostedMuon1;
	delete boostedMuon2;
	continue; 
      }
      if (kMuonRapRange      && (fMaxMuonRap<rap2      || rap2<fMinMuonRap)){
	delete boostedMuon1;
	delete boostedMuon2;
	continue;
      }
      if (kMuonEtaRange      && (fMaxMuonEta<eta2      || eta2<fMinMuonEta)){
	delete boostedMuon1;
	delete boostedMuon2;
	continue;
      }

      //printf("muon1 Kinematics: momentum = %.3f [GeV/c]   eta = %.3f\n",boostedMuon1->P(),boostedMuon1->PseudoRapidity());
      //printf("muon2 Kinematics: momentum = %.3f [GeV/c]   eta = %.3f\n",boostedMuon2->P(),boostedMuon2->PseudoRapidity());

      //Select charge
      if (gRandom->Rndm() < 0.5){
	pdgCode1 =  13;
	pdgCode2 = -13;
      }
      else{
	pdgCode1 = -13;
	pdgCode2 =  13;
      }

      //Store the muon pairs
      PushTrack(1, -1, Int_t(pdgCode1),
		px1,py1,pz1,energy1,
		origin[0],origin[1],origin[2],Double_t(time),
		polar[0],polar[1],polar[2],
		kPPrimary, nt, 1., 1);

      PushTrack(1, -1, Int_t(pdgCode2),
		px2, py2, pz2, energy2,
		origin[0], origin[1], origin[2], Double_t(time),
		polar[0], polar[1], polar[2],
		kPPrimary, nt, 1., 1);

      delete boostedMuon1;
      delete boostedMuon2;
      
      break;

    }
    
    
    break;
  }

  AliGenEventHeader* header = new AliGenEventHeader("ParamMuons");
  header->SetPrimaryVertex(fVertex);
  header->SetNProduced(nPartGenerated);
  header->SetInteractionTime(fTime);

  // Passes header either to the container or to gAlice
  if (fContainer) {
    fContainer->AddHeader(header);
  } else {
    gAlice->SetGenEventHeader(header);
  }
}

//====================================================================================================================================================

void AliGenDimuon::Init() {

  // Initialisation, check consistency of selected ranges

  if (TestBit(kPtRange) && TestBit(kMomentumRange))
    printf("You should not set the momentum range and the pt range at the same time!\n");
  if ((!TestBit(kPtRange)) && (!TestBit(kMomentumRange)))
    printf("You should set either the momentum or the pt range!\n");
  if ((TestBit(kYRange) && TestBit(kThetaRange)) || (TestBit(kYRange) && TestBit(kEtaRange)) || (TestBit(kEtaRange) && TestBit(kThetaRange)))
    printf("You should only set the range of one of these variables: y, eta or theta\n");
  if ((!TestBit(kYRange)) && (!TestBit(kEtaRange)) && (!TestBit(kThetaRange)))
    printf("You should set the range of one of these variables: y, eta or theta\n");
  
  printf("Parent PdgCode (%d)\n",fPdgCode);
  if(TestBit(kYRange)) 
    printf("Parent particle rapigity range %.3f<y<%.3f \n",fYMin,fYMax);
  if(TestBit(kPtRange)) 
    printf("Parent particle pT range %.3f<pT<%.3f \n",fPtMin,fPtMax);
  if(TestBit(kThetaRange)) 
    printf("Parent particle theta range %.3f<theta<%.3f \n",fThetaMin,fThetaMax);


  TFile* inFile = TFile::Open("./include/inputDimuonParam.root");
  
  if( fPdgCode == 113 ){
    fHistParentKine = (TH2F*)inFile->Get("fHistRhoRapPt");
  }
  else if( fPdgCode == 223 ){
    fHistParentKine = (TH2F*)inFile->Get("fHistOmegaRapPt");
  }
  else if( fPdgCode == 333 ){
    fHistParentKine = (TH2F*)inFile->Get("fHistPhiRapPt");
  }
  else if( fPdgCode == 443 ){
    fHistParentKine = (TH2F*)inFile->Get("fHistJpsiRapPt");
  }
  else if( fPdgCode ==  100443){
    fHistParentKine = (TH2F*)inFile->Get("fHistPsi2SRapPt");
  }
  else if( fPdgCode == 533 ){
    fHistParentKine = (TH2F*)inFile->Get("fHistUpsilon1SRapPt");
  }
  else if( fPdgCode == 100553 ){
    fHistParentKine = (TH2F*)inFile->Get("fHistUpsilon2SRapPt");
  }
  else if( fPdgCode == 200553 ){
    fHistParentKine = (TH2F*)inFile->Get("fHistUpsilon3SRapPt");
  }
  else{
    fHistParentKine = (TH2F*)inFile->Get("fHistPhiRapPt");
    printf("You should set the particle (113,223,333,433,100433,533,100533,200533). Now phi-meson (333) kinematic distribution is set, instead\n");
  }

  AliPDG::AddParticlesToPdgDataBase();
}

//====================================================================================================================================================

void AliGenDimuon::SetPtShape(){

  Double_t Ae = 187.;
  Double_t Te = 0.39;
  Double_t m0 = 0.135;
  Double_t A  = 1526.;
  Double_t T  = 0.29;
  Double_t n  = 2.75;

  fPt = new TF1("fPt","[0]*exp(-(sqrt(x*x + [1]*[1])-[1])/[2]) + [3]*pow(1+x*x/([4]*[4]*[5]),-[5])",0,10);
  fPt->SetParameters(Ae,Te,m0,A,T,n);
}

//====================================================================================================================================================
void AliGenDimuon::SetRapidityShape(){

  Double_t p0 = 1.87732;
  Double_t p1 = 0.00658212;
  Double_t p2 = -0.0988071;
  Double_t p3 = -0.000452746;
  Double_t p4 = 0.00269782;

  fRap = new TF1("fRap","[0]*(1 + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x)",-10,10);
  fRap->SetParameters(p0,p1,p2,p3,p4);
}

Double_t AliGenDimuon::EtaToTheta(Double_t arg) {
  return (180. / TMath::Pi()) * 2. * atan(exp(-arg));
}

/*
//====================================================================================================================================================
void AliGenDimuon::SetDimuonKinematics(){
  
  //TFile* input = 

  fRap = new TF1("fRap","[0]*(1 + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x)",-10,10);
  fRap->SetParameters(p0,p1,p2,p3,p4);
}
*/
