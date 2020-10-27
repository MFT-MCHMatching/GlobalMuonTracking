#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "SimulationDataFormat/MCTrack.h"

#endif

void Read_Kine_O2() {

  // DataFormats/simulation/include/SimulationDataFormat/MCTrack.h
  using o2::MCTrack;

  TFile fileK("o2sim_Kine.root");
  TTree* kineTree = (TTree*)fileK.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);

  int nEvents = kineTree->GetEntries();
  printf("Number of events %d \n", nEvents);
  
  for (int iev = 0; iev < nEvents; ++iev) {
    kineTree->GetEntry(iev);
    int nMCTracks = mcTrkVec.size();
    printf("Event %d has %d MC tracks\n", iev, nMCTracks);
    int nt = 0;
    for (auto mcTrack : mcTrkVec) {
      if (mcTrack.isSecondary()) continue;
      printf("MCTrack ID %4d   PDG %4d   name %s   isSec %d   E %7.3f \n",
	     nt,
	     mcTrack.GetPdgCode(),
	     TDatabasePDG::Instance()->GetParticle(mcTrack.GetPdgCode())->GetName(),
	     mcTrack.isSecondary(),
	     mcTrack.GetEnergy());
      ++nt;
    }
  }
  
  fileK.Close();
  
}
