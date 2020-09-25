#include <iostream>
// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONVDigit.h"
#include "AliESDMuonTrack.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMpDEIterator.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpPad.h"
#include "AliMpDetElement.h"
#include "AliMpCathodType.h"
#include "TFile.h"
#include "TTree.h"
#include "AliESDEvent.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliRunLoader.h"
#include "AliStack.h"



using namespace std;

#include "include/tempMCHTrackGetter.h"

#pragma link C++ class tempMCHTrack+;
#pragma link C++ class std::vector<tempMCHTrack>+;

std::vector<tempMCHTrack> tempMCHTracks;


Int_t fNMCH;

TTree* treeMCH = new TTree("treeMCH","treeMCH");


void ConvertMCHESDTracks(string in_dir=""){

  string out_dir = in_dir;
  if(!TFile::Open(Form("%s/AliESDs.root",in_dir.c_str()))) return;

  TString ocdbPath = "$ALIROOT_OCDB_ROOT/OCDB";
  TString fAlignOCDBpath     = "";
  TString fRecoParamOCDBpath = "";

  AliCDBManager *cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  if (!fAlignOCDBpath.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fAlignOCDBpath.Data());
  if (!fRecoParamOCDBpath.IsNull()) cdbm->SetSpecificStorage("MUON/Calib/RecoParam",fRecoParamOCDBpath.Data());
  cdbm->SetRun(255173);



  // load geometry for track extrapolation to vertex and for checking hits are under pads in reconstructible tracks
  if (!AliGeomManager::GetGeometry()) {
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
  }

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
      std::cout << " " << iTr << " " ;

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
