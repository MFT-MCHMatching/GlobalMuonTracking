void sim(Int_t nev = 4) {

  // libraries required by geant321
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6_4_25");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libgeant321");

  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALIDPG_ROOT/include");
  gROOT->LoadMacro("./AliGenMimicPbPb.cxx++g");
  gROOT->LoadMacro("./AliGenDimuon.cxx++g");

  AliSimulation simulator;
  simulator.SetMakeSDigits("MUON");
  // simulator.SetMakeDigitsFromHits("ITS");
  simulator.SetWriteRawData("MUON", "raw.root", kTRUE);

  simulator.SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
                               Form("local://%s", gSystem->pwd()));

  simulator.SetRunHLT(""); // In case we do not have ancored production

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
