#include <iostream>
// MUON includes
#include "TFile.h"
#include "TTree.h"
#include "AliESDEvent.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "include/tempMCHTrackGetter.h"


using namespace std;


#pragma link C++ class tempMCHTrack+;
#pragma link C++ class std::vector<tempMCHTrack>+;

std::vector<tempMCHTrack> tempMCHTracks;


Int_t fNMCH;

TTree* treeMCH = new TTree("treeMCH","treeMCH");


void ConvertMCHESDTracks(string in_dir=""){
  //gSystem->Load("libpythia6_4_25");

  string out_dir = in_dir;
  if(!TFile::Open(Form("%s/AliESDs.root",in_dir.c_str()))) return;

  TFile* inFile = new TFile(Form("%s/AliESDs.root",in_dir.c_str()),"read");

  TTree *esdTree = (TTree*)inFile->Get("esdTree");

  AliESDEvent *ev = new AliESDEvent();
  ev->ReadFromTree(esdTree);

  AliRunLoader* rl = AliRunLoader::Open(Form("%s/galice.root",in_dir.c_str()));
  rl->LoadgAlice();
  //gAlice = rl->GetAliRun();
  rl->LoadKinematics();
  rl->LoadHeader();

  treeMCH->Branch("tempMCHTracks",&tempMCHTracks);

  TFile* output  = new TFile(Form("%s/tempMCHTracks.root",out_dir.c_str()),"recreate");

  Long64_t nEvents = rl->GetNumberOfEvents();

  for (Int_t iEv = 0; iEv < nEvents; iEv++) {

    std::cout << "#Evt: " << iEv << std::endl;
    fNMCH=0;

    tempMCHTracks.clear();

    esdTree->GetEvent(iEv);
    rl->GetEvent(iEv);

    Int_t nTracks = ev->GetNumberOfMuonTracks();

    AliStack * stack = rl->Stack();
    //Int_t nPart = stack->GetNtrack();

    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliESDMuonTrack* esdTrack = (AliESDMuonTrack*)ev->GetMuonTrack(iTr);

      if (!esdTrack->ContainTrackerData() || esdTrack->GetLabel()<0) continue;

      TParticle * part = stack->Particle(esdTrack->GetLabel());

      if(esdTrack->GetMatchTrigger()<1)                                      continue;
      //if(esdTrack->GetChi2() / esdTrack->GetNDF()>5.0)                             continue;
      if(-3.6>esdTrack->Eta() && esdTrack->Eta()>-2.1)                             continue;
      if(17.6>esdTrack->GetRAtAbsorberEnd() && esdTrack->GetRAtAbsorberEnd()>89.5) continue;

      tempMCHTrackGetter t(*esdTrack);
      tempMCHTrack thisTrack;
      t.update(thisTrack);
      thisTrack.fiEv = iEv;
      tempMCHTracks.push_back(thisTrack);

      ++fNMCH;
    }

    treeMCH->Fill();

  }

  output->WriteTObject(treeMCH);
  std::cout << std::endl << output->GetName() << std::endl;
}
