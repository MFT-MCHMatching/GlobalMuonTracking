void finalizeGlobalFwdAssessment()
{
    o2::globaltracking::GloFwdAssessment analyser(true);
    analyser.loadHistos(); // loads GlobalForwardAssessment.root produced by the DPL workflow
    analyser.finalizeAnalysis();
    TFile *fout = new TFile("GlobalForwardAssessmentFinalized.root", "RECREATE");
    TObjArray objarOut;
    analyser.getHistos(objarOut);
    objarOut.Write();
    fout->Close();
}