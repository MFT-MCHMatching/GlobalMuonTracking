#include "GlobalTracking/MatchGlobalFwdAssessment.h"

// Dictionaries for MatchGlobalFwdAssessment.h are needed to have this fully working on a macro. 
// In the meantime one can run the following code on a root shell to finalize assessment

o2::globaltracking::GloFwdAssessment analyser(true);
analyser.loadHistos(); // loads GlobalForwardAssessment.root produced by the DPL workflow
analyser.finalizeAnalysis();
TFile* fout = new TFile("GlobalForwardAssessmentFinalized.root", "RECREATE");
TObjArray objarOut;
analyser.getHistos(objarOut);
objarOut.Write();
fout->Close();