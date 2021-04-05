#include "AliGenParam.h"
#include "AliGenerator.h"
#include "TRandom.h"
#include <Riostream.h>

enum PprTrigConf_t { kDefaultPPTrig, kDefaultPbPbTrig };

const char *pprTrigConfName[] = {"p-p", "Pb-Pb"};

Float_t EtaToTheta(Float_t arg);
Int_t IpMuon(TRandom *ran);
Double_t PtMuon(const Double_t *px, const Double_t * /*dummy*/);
Double_t YMuon(const Double_t *py, const Double_t * /*dummy*/);
Double_t V2Muon(const Double_t * /*dummy*/, const Double_t * /*dummy*/);

Int_t IpPion(TRandom *ran);
Double_t PtPion(const Double_t *px, const Double_t * /*dummy*/);
Double_t YPion(const Double_t *py, const Double_t * /*dummy*/);
Double_t V2Pion(const Double_t * /*dummy*/, const Double_t * /*dummy*/);

static AliMagF::BeamType_t beamType = AliMagF::kBeamTypepp;
static Double_t beamEnergy = 7000;           //.*82./208;
static PprTrigConf_t strig = kDefaultPPTrig; // default PP trigger configuration

void Config() {
  // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
  // Theta range given through pseudorapidity limits 22/6/2001

  // Set Random Number seed
  gRandom->SetSeed(std::atoi(std::getenv("SEED"))); // Set 0 to use the currecnt time
  // The libraries required by Geant3 are loaded in sim.C

  new TGeant3TGeo("C++ Interface to Geant3");

  AliRunLoader *rl = 0x0;

  rl = AliRunLoader::Open("galice.root", AliConfig::GetDefaultEventFolderName(),
                          "recreate");
  if (rl == 0x0) {
    gAlice->Fatal("Config.C", "Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(10000);
  gAlice->SetRunLoader(rl);

  // Set the trigger configuration
  AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
  cout << "Trigger configuration is set to  " << pprTrigConfName[strig] << endl;

  //
  // Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();

  decayer->SetForceDecay(kAll);
  decayer->Init();

  TVirtualMC *vmc = TVirtualMC::GetMC();

  vmc->SetExternalDecayer(decayer);
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV

  vmc->SetProcess("DCAY", 1);
  vmc->SetProcess("PAIR", 1);
  vmc->SetProcess("COMP", 1);
  vmc->SetProcess("PHOT", 1);
  vmc->SetProcess("PFIS", 0);
  vmc->SetProcess("DRAY", 0);
  vmc->SetProcess("ANNI", 1);
  vmc->SetProcess("BREM", 1);
  vmc->SetProcess("MUNU", 1);
  vmc->SetProcess("CKOV", 1);
  vmc->SetProcess("HADR", 1);
  vmc->SetProcess("LOSS", 2);
  vmc->SetProcess("MULS", 1);
  vmc->SetProcess("RAYL", 1);

  Float_t cut = 1.e-3; // 1MeV cut by default
  Float_t tofmax = 1.e10;

  vmc->SetCut("CUTGAM", cut);
  vmc->SetCut("CUTELE", cut);
  vmc->SetCut("CUTNEU", cut);
  vmc->SetCut("CUTHAD", cut);
  vmc->SetCut("CUTMUO", cut);
  vmc->SetCut("BCUTE", cut);
  vmc->SetCut("BCUTM", cut);
  vmc->SetCut("DCUTE", cut);
  vmc->SetCut("DCUTM", cut);
  vmc->SetCut("PPCUTM", cut);
  vmc->SetCut("TOFMAX", tofmax);

  // Special generation for Valgrind tests
  // Each detector is fired by few particles selected
  // to cover specific cases

  AliGenCocktail *gener = new AliGenCocktail();
  gener->SetEnergyCMS(beamEnergy); // Needed by ZDC
  gener->SetPhiRange(0, 360);
  // Set pseudorapidity range from -2.2 to -3.7
  Float_t thmin = EtaToTheta(-2.2); // theta min. <---> eta max
  Float_t thmax = EtaToTheta(-3.7); // theta max. <---> eta min
  gener->SetThetaRange(thmin, thmax);
  gener->SetOrigin(0, 0, 0); // vertex position
  gener->SetSigma(0, 0, 0);  // Sigma in (X,Y,Z) (cm) on IP position

  std::string MCHgen;
  if (gSystem->Getenv("MCHGENERATOR")) {
    MCHgen = gSystem->Getenv("MCHGENERATOR");
    std::cout << " Using MCHGENERATOR: " << MCHgen << std::endl;
  } else {
    MCHgen = "gun0_100GeV";
    std::cout << " No MCHGENERATOR defined. Using default: " << MCHgen
              << std::endl;
  }

  ofstream genMatcherLog("MatcherGenConfig.txt");
  genMatcherLog << MCHgen;
  // Generators
  if (MCHgen.find("gun0_100GeV") < MCHgen.length()) {
    std::cout << " This is gun0_100GeV Generator! " << std::endl;

    Int_t nPions = 10;
    if (gSystem->Getenv("NPIONS")) {
      nPions = atoi(gSystem->Getenv("NPIONS"));
      nPions = 2 * (nPions / 2);
    }

    Int_t nMuons = 2;
    if (gSystem->Getenv("NMUONS")) {
      nMuons = atoi(gSystem->Getenv("NMUONS"));
      nMuons = 2 * (nMuons / 2);
    }

    std::cout << " nPions =  " << nPions << std::endl;
    std::cout << " nMuons =  " << nMuons << std::endl;

    genMatcherLog << "_" << nPions << "pi_" << nMuons << "mu_";

    // Pions
    AliGenBox *gPPions = new AliGenBox(nPions / 2);
    gPPions->SetMomentumRange(0.1, 100.1);
    gPPions->SetPhiRange(0., 360.);
    gPPions->SetThetaRange(thmin, thmax);
    gPPions->SetPart(kPiPlus); // Positive pions
    gener->AddGenerator(gPPions, "POS PIONS", 1);

    AliGenBox *gNPions = new AliGenBox(nPions / 2);
    gNPions->SetMomentumRange(0.1, 100.1);
    gNPions->SetPhiRange(0., 360.);
    gNPions->SetThetaRange(thmin, thmax);
    gNPions->SetPart(kPiMinus); // Positive pions
    gener->AddGenerator(gNPions, "NEG PIONS", 1);

    // MUONS
    AliGenBox *gmuon1 = new AliGenBox(nMuons / 2);
    gmuon1->SetMomentumRange(0.1, 100);
    gmuon1->SetPhiRange(0., 360.);
    gmuon1->SetThetaRange(thmin, thmax);
    gmuon1->SetPart(kMuonMinus); // Negative muons
    gener->AddGenerator(gmuon1, "GENBOX MUON1", 1);

    AliGenBox *gmuon2 = new AliGenBox(nMuons / 2);
    gmuon2->SetMomentumRange(0.1, 100);
    gmuon2->SetPhiRange(0., 360.);
    gmuon2->SetThetaRange(thmin, thmax);
    gmuon2->SetPart(kMuonPlus); // Positive muons
    gener->AddGenerator(gmuon2, "GENBOX MUON2", 1);
  }
  //
  // Muon Box Gun
  if (MCHgen.find("MuBoxGun") < MCHgen.length()) {
    std::cout << " This is Mu_gun0_100GeV Generator! " << std::endl;

    Int_t nMuons = 2;
    if (gSystem->Getenv("NMUONS")) {
      nMuons = atoi(gSystem->Getenv("NMUONS"));
      nMuons = 2 * (nMuons / 2);
    }

    std::cout << " nMuons =  " << nMuons << std::endl;

    genMatcherLog << "_" << nMuons << "mu_";

    // MUONS
    AliGenBox *gmuon1 = new AliGenBox(nMuons / 2);
    gmuon1->SetMomentumRange(0.1, 100);
    gmuon1->SetPhiRange(0., 360.);
    gmuon1->SetThetaRange(thmin, thmax);
    gmuon1->SetPart(kMuonMinus); // Negative muons
    gener->AddGenerator(gmuon1, "GENBOX NMUONS", 1);

    AliGenBox *gmuon2 = new AliGenBox(nMuons / 2);
    gmuon2->SetMomentumRange(0.1, 100);
    gmuon2->SetPhiRange(0., 360.);
    gmuon2->SetThetaRange(thmin, thmax);
    gmuon2->SetPart(kMuonPlus); // Positive muons
    gener->AddGenerator(gmuon2, "GENBOX PMUONS", 1);
  }

  if (MCHgen.find("PiMuParam") < MCHgen.length()) {
    std::cout << " This is PiMuParam generator! " << std::endl;

    Int_t nPions = 20;
    if (gSystem->Getenv("NPIONS")) {
      nPions = atoi(gSystem->Getenv("NPIONS"));
    }
    std::cout << " nPions = " << nPions << std::endl;
    Int_t nMuons = 2;
    if (gSystem->Getenv("NMUONS")) {
      nMuons = atoi(gSystem->Getenv("NMUONS"));
    }
    std::cout << " nMuons = " << nMuons << std::endl;
    std::cout << " nPions = " << nPions << std::endl;

    genMatcherLog << "_" << nPions << "pi_" << nMuons << "mu_";

    // Pions
    AliGenParam *gPions =
        new AliGenParam(nPions, -1, PtPion, YPion, V2Pion, IpPion);
    gPions->SetMomentumRange(0., 1.e6);
    gPions->SetPtRange(0.1, 999.);
    gPions->SetYRange(-4.3, -2.3);
    gPions->SetPhiRange(0., 360.);
    gPions->SetForceDecay(
        kNoDecay); // Why is this needed to generate MFT tracks in O2?
    gPions->SetTrackingFlag(1);
    gener->AddGenerator(gPions, "PIONS", 1);

    // Muons
    AliGenParam *gMuons =
        new AliGenParam(nMuons, -1, PtMuon, YMuon, V2Muon, IpMuon);
    gMuons->SetMomentumRange(4., 1.e6);
    gMuons->SetPtRange(0.0, 999.);
    gMuons->SetYRange(-4.3, -2.3);
    gMuons->SetPhiRange(0., 360.);
    gMuons->SetForceDecay(kNoDecay);
    gMuons->SetTrackingFlag(1);
    gener->AddGenerator(gMuons, "MUONS", 1);
  }

  //
  if (MCHgen.find("PiParam") < MCHgen.length()) {
    std::cout << " This is PiParam generator! " << std::endl;

    Int_t nPions = 20;
    if (gSystem->Getenv("NPIONS")) {
      nPions = atoi(gSystem->Getenv("NPIONS"));
      std::cout << " Defined nPions =  " << nPions << std::endl;
    }
    std::cout << " nPions = " << nPions << std::endl;
    genMatcherLog << "_" << nPions << "pi";

    // Pions
    AliGenParam *gPions =
        new AliGenParam(nPions, -1, PtPion, YPion, V2Pion, IpPion);
    gPions->SetMomentumRange(0., 1.e6);
    gPions->SetPtRange(0.1, 999.);
    gPions->SetYRange(-4.3, -2.3);
    gPions->SetPhiRange(0., 360.);
    gPions->SetForceDecay(
        kNoDecay); // Why is this needed to generate MFT tracks in O2?
    gPions->SetTrackingFlag(1);
    gener->AddGenerator(gPions, "PIONS", 1);
  }

  if (MCHgen.find("JPsiParam") < MCHgen.length()) {
    std::cout << " This is JPsiParam generator! " << std::endl;

    Int_t nJPsi;
    if (gSystem->Getenv("NJPSI")) {
      nJPsi = atoi(gSystem->Getenv("NJPSI"));
      std::cout << " Defined nJPsi =  " << nJPsi << std::endl;

    } else {
      std::cout << " Default nJPsi = 1 " << std::endl;
      nJPsi = 1;
    }

    genMatcherLog << "_" << nJPsi << "jpsi";

    AliGenParam *gJPsi = new AliGenParam(nJPsi, AliGenMUONlib::kJpsi);
    gJPsi->SetMomentumRange(0, 999);
    gJPsi->SetPtRange(0, 100.);
    gJPsi->SetPhiRange(0., 360.);
    gJPsi->SetCutOnChild(1);
    gJPsi->SetChildPhiRange(0., 360.);
    gJPsi->SetChildThetaRange(171.0, 178.0);
    gJPsi->SetOrigin(0, 0, 0);
    gJPsi->SetForceDecay(kDiMuon);
    gJPsi->SetTrackingFlag(1);
    gener->AddGenerator(gJPsi, "JPSI", 1);
  }

  if (MCHgen.find("hijing") < MCHgen.length()) { // Hijing generator
    std::cout << " This is hijing generator! " << std::endl;
    /*
    AliGenHijing *gener = new AliGenHijing(-1);
    // centre of mass energy
    gener->SetEnergyCMS(5500.);
    // reference frame
    gener->SetReferenceFrame("CMS");
    // projectile
    gener->SetProjectile("A", 208, 82);
    gener->SetTarget("A", 208, 82);
    // tell hijing to keep the full parent child chain
    gener->KeepFullEvent();
    // enable jet quenching
    gener->SetJetQuenching(1);
    // enable shadowing
    gener->SetShadowing(1);
    // neutral pion and heavy particle decays switched off
    gener->SetDecaysOff(1);
    // Don't track spectators
    gener->SetSpectators(0);
    // kinematic selection
    gener->SetSelectAll(0);
    // impact parameter range
    gener->SetImpactParameterRange(0., 5.); // 0. - 5. fm corresponds to ~10% most central
    gener->Init();
    */

    AliGenMimicPbPb *mimic_pbpb = new AliGenMimicPbPb(1);
    mimic_pbpb->Init();
    gener->AddGenerator(mimic_pbpb, "MIMIC_PBPB", 1);
    
    
  }
  if (MCHgen.find("muoncocktail") < MCHgen.length()) { // Muon cocktail for PbPb
    std::cout << " This is muoncocktail generator! " << std::endl;
    AliGenMUONCocktail *generMUONCocktail = new AliGenMUONCocktail();
    generMUONCocktail->SetPtRange(1., 100.);  // Transverse momentum range
    generMUONCocktail->SetPhiRange(0., 360.); // Azimuthal angle range
    generMUONCocktail->SetYRange(-4.0, -2.5);
    generMUONCocktail->SetMuonPtCut(0.5);
    generMUONCocktail->SetMuonThetaCut(171., 178.);
    generMUONCocktail->SetMuonMultiplicity(2);
    generMUONCocktail->SetImpactParameterRange(
        0., 5.); // 10% most centra PbPb collisions
    generMUONCocktail->SetVertexSmear(kPerTrack);
    generMUONCocktail->SetOrigin(0, 0, 0); // Vertex position
    generMUONCocktail->SetSigma(0, 0,
                                0.0); // Sigma in (X,Y,Z) (cm) on IP position
    gener->AddGenerator(generMUONCocktail, "MUONCocktail", 1);
  }

  Int_t nEvts;
  if (gSystem->Getenv("NEV")) {
    nEvts = atoi(gSystem->Getenv("NEV"));
  }
  genMatcherLog << "_" << nEvts << "evts" << std::endl;
  genMatcherLog.close();

  gener->Init();

  //
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //
  // gener->SetVertexSmear(perTrack);
  // Field (L3 0.5 T)
  TGeoGlobalMagField::Instance()->SetField(new AliMagF(
      "Maps", "Maps", -1., -1., AliMagF::k5kG, beamType, beamEnergy));

  Int_t iABSO = 1;
  Int_t iDIPO = 1;
  Int_t iFMD = 0;
  Int_t iFRAME = 1;
  Int_t iHALL = 1;
  Int_t iITS = 0;
  Int_t iMAG = 1;
  Int_t iMUON = 1;
  Int_t iPHOS = 0;
  Int_t iPIPE = 1;
  Int_t iPMD = 0;
  Int_t iHMPID = 0;
  Int_t iSHIL = 1;
  Int_t iT0 = 0;
  Int_t iTOF = 0;
  Int_t iTPC = 0;
  Int_t iTRD = 0;
  Int_t iZDC = 0;
  Int_t iEMCAL = 0;
  Int_t iACORDE = 0;
  Int_t iVZERO = 0;
  rl->CdGAFile();
  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");

  if (iMAG) {
    //=================== MAG parameters ============================
    // --- Start with Magnet since detector layouts may be depending ---
    // --- on the selected Magnet dimensions ---
    AliMAG *MAG = new AliMAG("MAG", "Magnet");
  }

  if (iABSO) {
    //=================== ABSO parameters ============================
    AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
  }

  if (iDIPO) {
    //=================== DIPO parameters ============================

    AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 3");
  }

  if (iHALL) {
    //=================== HALL parameters ============================

    AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
  }

  if (iFRAME) {
    //=================== FRAME parameters ============================

    AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
    FRAME->SetHoles(1);
  }

  if (iSHIL) {
    //=================== SHIL parameters ============================

    AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 3");
  }

  if (iPIPE) {
    //=================== PIPE parameters ============================

    AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
  }

  if (iITS) {
    //=================== ITS parameters ============================

    AliITS *ITS = new AliITSv11("ITS", "ITS v11");
  }

  if (iTPC) {
    //============================ TPC parameters ===================
    AliTPC *TPC = new AliTPCv2("TPC", "Default");
  }

  if (iTOF) {
    //=================== TOF parameters ============================
    AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
  }

  if (iHMPID) {
    //=================== HMPID parameters ===========================
    AliHMPID *HMPID = new AliHMPIDv3("HMPID", "normal HMPID");
  }

  if (iZDC) {
    //=================== ZDC parameters ============================

    AliZDC *ZDC = new AliZDCv4("ZDC", "normal ZDC");
  }

  if (iTRD) {
    //=================== TRD parameters ============================

    AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
    AliTRDgeometry *geoTRD = TRD->GetGeometry();
    // Partial geometry: modules at 0,1,7,8,9,10,17
    // starting at 3h in positive direction
    geoTRD->SetSMstatus(2, 0);
    geoTRD->SetSMstatus(3, 0);
    geoTRD->SetSMstatus(4, 0);
    geoTRD->SetSMstatus(5, 0);
    geoTRD->SetSMstatus(6, 0);
    geoTRD->SetSMstatus(11, 0);
    geoTRD->SetSMstatus(12, 0);
    geoTRD->SetSMstatus(13, 0);
    geoTRD->SetSMstatus(14, 0);
    geoTRD->SetSMstatus(15, 0);
    geoTRD->SetSMstatus(16, 0);
  }

  if (iFMD) {
    //=================== FMD parameters ============================
    AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
  }

  if (iMUON) {
    //=================== MUON parameters ===========================
    // New MUONv1 version (geometry defined via builders)
    AliMUON *MUON = new AliMUONv1("MUON", "default");
  }
  //=================== PHOS parameters ===========================

  if (iPHOS) {
    AliPHOS *PHOS = new AliPHOSv1("PHOS", "Run1");
  }

  if (iPMD) {
    //=================== PMD parameters ============================
    AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
  }

  if (iT0) {
    //=================== T0 parameters ============================
    AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
  }

  if (iEMCAL) {
    //=================== EMCAL parameters ============================
    AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
  }

  if (iACORDE) {
    //=================== ACORDE parameters ============================
    AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
  }

  if (iVZERO) {
    //=================== VZERO parameters ============================
    AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
  }
}

Float_t EtaToTheta(Float_t arg) {
  return (180. / TMath::Pi()) * 2. * atan(exp(-arg));
}

Int_t IpMuon(TRandom *ran) {
  // muon composition

  if (ran->Rndm() < 0.5) {
    return 13;
  } else {
    return -13;
  }
}

//-------------------------------------------------------------------------
Double_t PtMuon(const Double_t *px, const Double_t * /*dummy*/) {
  // muon pT

  Double_t x = *px;

  Double_t p0 = 797.446;
  Double_t p1 = 0.830278;
  Double_t p2 = 0.632177;
  Double_t p3 = 10.2202;
  Double_t p4 = -0.000614809;
  Double_t p5 = -1.70993;

  return p0 * (1. / TMath::Power(p1 + TMath::Power(x, p2), p3) +
               p4 * TMath::Exp(p5 * x));
}

//-------------------------------------------------------------------------
Double_t YMuon(const Double_t *py, const Double_t * /*dummy*/) {
  // muon y

  Double_t y = *py;

  Double_t p0 = 1.87732;
  Double_t p1 = 0.00658212;
  Double_t p2 = -0.0988071;
  Double_t p3 = -0.000452746;
  Double_t p4 = 0.00269782;

  return p0 * (1. + p1 * y + p2 * y * y + p3 * y * y * y + p4 * y * y * y * y);
}

//-------------------------------------------------------------------------
Double_t V2Muon(const Double_t * /*dummy*/, const Double_t * /*dummy*/) {
  // muon v2
  return 0.;
}

//-------------------------------------------------------------------------
Int_t IpPion(TRandom *ran) {
  // Pion composition

  if (ran->Rndm() < 0.5) {
    return 211;
  } else {
    return -211;
  }
}

//-------------------------------------------------------------------------
Double_t PtPion(const Double_t *px, const Double_t * /*dummy*/) {
  // Pion pT

  Double_t x = *px;

  Double_t p0 = 797.446;
  Double_t p1 = 0.830278;
  Double_t p2 = 0.632177;
  Double_t p3 = 10.2202;
  Double_t p4 = -0.000614809;
  Double_t p5 = -1.70993;

  return p0 * (1. / TMath::Power(p1 + TMath::Power(x, p2), p3) +
               p4 * TMath::Exp(p5 * x));
}

//-------------------------------------------------------------------------
Double_t YPion(const Double_t *py, const Double_t * /*dummy*/) {
  // Pion y

  Double_t y = *py;

  Double_t p0 = 1.87732;
  Double_t p1 = 0.00658212;
  Double_t p2 = -0.0988071;
  Double_t p3 = -0.000452746;
  Double_t p4 = 0.00269782;

  return p0 * (1. + p1 * y + p2 * y * y + p3 * y * y * y + p4 * y * y * y * y);
}

//-------------------------------------------------------------------------
Double_t V2Pion(const Double_t * /*dummy*/, const Double_t * /*dummy*/) {
  // Pion v2
  return 0.;
}
