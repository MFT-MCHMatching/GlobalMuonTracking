void finalizeMFTAssessment()
{
  o2::mft::MFTAssessment analyser(true);
  analyser.loadHistos(); // loads MFFAssessment.root produced by the DPL workflow
  analyser.finalizeAnalysis();
  TFile* fout = new TFile("MFTAssessmentFinalized.root", "RECREATE");
  TObjArray objarOut;
  analyser.getHistos(objarOut); // Output objects
  objarOut.Write();
  fout->Close();
}