#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TParticle.h>

#endif

void Read_Kine_AR() {

  // Run loader
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");

  // gAlice
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();

  // Load kinematics and event header
  runLoader->LoadKinematics();
  runLoader->LoadHeader();

  // Particle properties
  TDatabasePDG * pdgdb = TDatabasePDG::Instance();

  Int_t nevMC = (Int_t)runLoader->GetNumberOfEvents();
  printf("Events: %d \n", nevMC);
  for (Int_t iev = 0; iev < nevMC; iev++) {
    // Get MC event
    runLoader->GetEvent(iev);
    // Particle stack
    AliStack * stack = runLoader->Stack();
    Int_t npart = stack->GetNtrack();
    printf("Events %d particles in stack %d \n", iev, npart);
    Int_t nt = 0;
    for (Int_t ipart = 0; ipart < npart; ipart++) {
      TParticle * part = stack->Particle(ipart);
      Int_t pdgcode = part->GetPdgCode();
      TParticlePDG * pdgpart = pdgdb->GetParticle(pdgcode);
      Int_t mumid = part->GetFirstMother();
      if (mumid < 0) {
	//pdgpart->Print();
	printf("MCTrack ID %4d   PDG %4d   name %s   isSec %d   E %7.3f \n",
	       nt,
	       part->GetPdgCode(),
	       part->GetName(),
	       (mumid < 0) ? 0 : 1,
	       part->Energy());
	++nt;
      }
    }

  }

  delete runLoader;

}

