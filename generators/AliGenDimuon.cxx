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

AliGenDimuon::AliGenDimuon():
  AliGenerator(), 
  fPt(0x0),
  fRap(0x0),
  fPdgCode(0),
  fBW(0x0){

  // Default constructor
    
}

//====================================================================================================================================================

AliGenDimuon::AliGenDimuon(Int_t nPart/*, Char_t *inputFile*/):
  AliGenerator(nPart),
  fPt(0x0),
  fRap(0x0),
  fPdgCode(0){

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
  Double_t rap=0.;
  Double_t mom=0.;
  Double_t energy=0;
  Double_t phi=0.;
  Double_t time=0.;
  Double_t theta = 0.;
    
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
  
  Double_t energy2 = 0;
  Double_t px2     = 0;
  Double_t py2     = 0;
  Double_t pz2     = 0;
  Double_t pt2     = 0;
  Double_t mom2    = 0;
  Double_t rap2    = 0;
  Double_t theta2  = 0;
 
  for (Int_t j=0; j<3; j++) origin[j] = fOrigin[j];
  time = fTimeOrigin;
  if (fVertexSmear==kPerEvent) {
    Vertex();
    for (Int_t j=0; j<3; j++) origin[j] = fVertex[j];
    time = fTime;
  }

  Int_t nPartGenerated = 2;
    
  Double_t m_muon       = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  Double_t mass_parent  = TDatabasePDG::Instance()->GetParticle(fPdgCode)->Mass();
  Double_t width_parent = TDatabasePDG::Instance()->GetParticle(fPdgCode)->Width();
  fBW->SetParameters(1,mass_parent,width_parent);

  mass = mass_parent;

  /*
  while(1){

    //mass = fBW->GetRandom();
    mass = mass_parent;
    
    if(mass<2*m_muon){
      continue;
    }
    else{
      break;
    }        
  }
  */
  
  Double_t daughter_m[2]={m_muon,m_muon};
  
  TLorentzVector rest_p(0,0,0,mass);

  TGenPhaseSpace ps_decay;
  ps_decay.SetDecay(rest_p,2,daughter_m);
  
  while (1){;
    
    ps_decay.Generate();
    
    pt  = fPt  -> GetRandom(0.,10);
    rap = fRap -> GetRandom(-3.6,-2.5);
    phi = gRandom->Uniform(0.,TMath::TwoPi());
    
    if (TestBit(kPtRange)       && (pt<fPtMin || pt>fPtMax))             continue;
    if (TestBit(kYRange)        && (rap<fYMin || rap>fYMax))             continue;
  
    TLorentzVector parent;
    parent.SetPtEtaPhiM(pt,rap,phi,mass);
    
    mom   = parent.P();
    theta = parent.Theta();
    
    TVector3 vec_boost =parent.BoostVector();
    
    TLorentzVector* muon1 = ps_decay.GetDecay(0);
    TLorentzVector* muon2 = ps_decay.GetDecay(1);

    muon1->Boost(vec_boost);
    muon2->Boost(vec_boost);
        
    energy1 = muon1->E();
    px1     = muon1->Px();
    py1     = muon1->Py();
    pz1     = muon1->Pz();
    pt1     = muon1->Pt();
    mom1    = muon1->P();
    rap1    = muon1->Rapidity();
    theta1  = muon1->Theta();

    //if (TestBit(kPtRange)       && (pt1<fPtMin || pt1>fPtMax))             continue;
    if (TestBit(kYRange)        && (rap1<fYMin || rap1>fYMax))             continue;
    //if (TestBit(kMomentumRange) && (mom1<fPMin || mom1>fPMax))             continue;    
    //if (TestBit(kThetaRange)    && (theta1<fThetaMin || theta1>fThetaMax)) continue;

    if(mom1<4.0) continue;
    
    energy2 = muon2->E();
    px2     = muon2->Px();
    py2     = muon2->Py();
    pz2     = muon2->Pz();
    pt2     = muon2->Pt();
    mom2    = muon2->P();
    rap2    = muon2->Rapidity();
    theta2  = muon2->Theta();

    //if (TestBit(kPtRange)       && (pt2<fPtMin || pt2>fPtMax))             continue;
    if (TestBit(kYRange)        && (rap2<fYMin || rap2>fYMax))             continue;
    //if (TestBit(kMomentumRange) && (mom2<fPMin || mom2>fPMax))             continue;    
    //if (TestBit(kThetaRange)    && (theta2<fThetaMin || theta2>fThetaMax)) continue;

    if(mom2<4.0) continue;
    
    if (gRandom->Rndm() < 0.5){
      pdgCode1 =  13;
      pdgCode2 = -13;
    }
    else{
      pdgCode1 = -13;
      pdgCode2 =  13;
    }

    PushTrack(1, -1, Int_t(pdgCode1),
              px1,py1,pz1,energy1,
              origin[0],origin[1],origin[2],Double_t(time),
              polar[0],polar[1],polar[2],
              kPPrimary, nt, 1., 1);
    
    PushTrack(1, -1, Int_t(pdgCode2),
              px2,py2,pz2,energy2,
              origin[0],origin[1],origin[2],Double_t(time),
              polar[0],polar[1],polar[2],
              kPPrimary, nt, 1., 1);    

    //cout<<pt<<"    "<<mom1<<"    "<<mom2<<endl;
    
    TLorentzVector muon12 = *muon1 + *muon2;
    //cout<<muon12.Pt()<<"   "<<muon12.Eta()<<"   "<<muon12.M()<<endl;
    
    /*
    Double_t m1  = muon1->M();
    Double_t px1 = muon1->Px();
    Double_t py1 = muon1->Py();
    Double_t pz1 = muon1->Pz();
    Double_t p1  = sqrt(px1*px1 + py1*py1 + pz1*pz1);    
    Double_t e1  = sqrt(m1*m1 + p1*p1);

    Double_t m2  = muon2->M();
    Double_t px2 = muon2->Px();
    Double_t py2 = muon2->Py();
    Double_t pz2 = muon2->Pz();
    Double_t p2  = sqrt(px2*px2 + py2*py2 + pz2*pz2);    
    Double_t e2  = sqrt(m2*m2 + p2*p2);
    */
    
    break;
  }

  AliGenEventHeader* header = new AliGenEventHeader("ParamMuons");
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

void AliGenDimuon::Init() {

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
  ///*
  if (TestBit(kPtRange) && TestBit(kMomentumRange)) 
    printf("You should not set the momentum range and the pt range at the same time!\n");
  if ((!TestBit(kPtRange)) && (!TestBit(kMomentumRange))) 
    printf("You should set either the momentum or the pt range!\n");
  if ((TestBit(kYRange) && TestBit(kThetaRange)) || (TestBit(kYRange) && TestBit(kEtaRange)) || (TestBit(kEtaRange) && TestBit(kThetaRange)))
    printf("You should only set the range of one of these variables: y, eta or theta\n");
  if ((!TestBit(kYRange)) && (!TestBit(kEtaRange)) && (!TestBit(kThetaRange)))
    printf("You should set the range of one of these variables: y, eta or theta\n");
  //*/
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
