#if !defined(__CLING__) || defined(__ROOTCLING__)

#ifdef __MAKECINT__
//#pragma link C++ class GlobalMuonTrack + ;
//#pragma link C++ class std::vector < GlobalMuonTrack> + ;
//#pragma link C++ class MatchingHelper + ;
#endif

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "ReconstructionDataFormats/GlobalFwdTrack.h"
#include "DataFormatsMFT/TrackMFT.h"

#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TProfile.h>
#include <TStyle.h>

#endif

#include "macrohelpers/HistosHelpers.C"
#include "macrohelpers/MagField.C"

using o2::MCTrackT;
using GlobalMuonTrack = o2::dataformats::GlobalFwdTrack;
using eventFoundTracks = std::vector<bool>;
using std::vector;

bool DEBUG_VERBOSE = false;
bool EXPORT_HISTOS_IMAGES = false;

//_________________________________________________________________________________________________
int GlobalMuonChecks(const std::string trkFile = "globalfwdtracks.root",
                     const std::string o2sim_KineFile = "o2sim_Kine.root",
                     const std::string sig_KineFile = "sgn_Kine.root")
{

  if (gSystem->Getenv("VERBOSEMATCHING")) {
    std::cout << " Vebose checking enabled." << std::endl;
    DEBUG_VERBOSE = true;
  }

  // Histos parameters
  Double_t pMin = 0.0;
  Double_t pMax = 100.0;
  Double_t deltaetaMin = -.1;
  Double_t deltaetaMax = +.1;
  Double_t etaMin = -3.5;
  Double_t etaMax = -2.4;
  Double_t deltaphiMin = -.01;
  Double_t deltaphiMax = .01;
  Double_t deltatanlMin = -0.2;
  Double_t deltatanlMax = 0.2;

  /*
  // histos
  // gROOT->SetStyle("Bold");
  gStyle->SetOptStat("emr");
  gStyle->SetStatW(.28);
  gStyle->SetStatH(.26);
  gStyle->SetPalette(1, 0);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  // gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(2);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.06, "xyz");
  gStyle->SetLabelOffset(0.01, "y");
  // gStyle->SetLabelColor(kBlue,"xy");
  gStyle->SetTitleSize(0.06, "xyz");
  gStyle->SetTitleSize(0.08, "o");
  gStyle->SetTitleOffset(0.95, "Y");
  gStyle->SetTitleFillColor(10);
  // gStyle->SetTitleTextColor(kNlacBlue);
  gStyle->SetStatColor(10);
*/
  enum TH2HistosCodes {
    kGMTrackDeltaXYVertex,
    kGMTrackQPtRec_MC,
    kGMTrackPtResolution,
    kGMTrackInvPtResolution,
    kMCTracksEtaZ
  };

  std::map<int, const char*> TH2Names{
    {kGMTrackDeltaXYVertex, "Global Muon Tracks Vertex at Z = 0"},
    {kGMTrackQPtRec_MC, "GM Track QPt FITxMC"},
    {kGMTrackPtResolution, "GM Track Pt Resolution"},
    {kGMTrackInvPtResolution, "GM Track InvPt Resolution"},
    {kMCTracksEtaZ, "MCTracks_eta_z"}};

  std::map<int, const char*> TH2Titles{
    {kGMTrackDeltaXYVertex, "Global Muon Tracks at Z_vertex"},
    {kGMTrackQPtRec_MC, "q*Pt: Reconstructed vs MC"},
    {kGMTrackPtResolution, "Pt Resolution"},
    {kGMTrackInvPtResolution, "InvPt Resolution"},
    {kMCTracksEtaZ, "MC Tracks: Pseudorapidity vs zVertex"}};

  std::map<int, std::array<double, 6>> TH2Binning{
    {kGMTrackDeltaXYVertex, {100, -.5, .5, 100, -.5, .5}},
    {kGMTrackQPtRec_MC, {40, -20, 20, 40, -20, 20}},
    {kGMTrackPtResolution, {40, 0, 20, 100, 0, 5}},
    {kGMTrackInvPtResolution, {14, 0, 7, 300, -2, 2}},
    {kMCTracksEtaZ, {31, -15, 16, 25, etaMin, etaMax}}};

  std::map<int, const char*> TH2XaxisTitles{
    {kGMTrackDeltaXYVertex, "\\Delta x ~[mm]"},
    {kGMTrackQPtRec_MC, "(q.pt)_{MC} [GeV]"},
    {kGMTrackPtResolution, "pt_{MC} [GeV]"},
    {kGMTrackInvPtResolution, "pt_{MC} [GeV]"},
    {kMCTracksEtaZ, "Vertex PosZ [cm]"}};

  std::map<int, const char*> TH2YaxisTitles{
    {kGMTrackDeltaXYVertex, "\\Delta y ~[mm]"},
    {kGMTrackQPtRec_MC, "(q.pt)_{fit} [GeV]"},
    {kGMTrackPtResolution, "pt_{fit} / pt_{MC}"},
    {kGMTrackInvPtResolution, "(1/(p_t)_{fit} - 1/(p_t)_{MC})*(p_t)_{MC}"},
    {kMCTracksEtaZ, "\\eta"}};

  enum TH1HistosCodes {
    kGMTrackPullX,
    kGMTrackPullY,
    kGMTrackPullPhi,
    kGMTrackPullTanl,
    kGMTrackPullInvQPt,
    kGMTrackMatchingPullX,
    kGMTrackMatchingPullY,
    kGMTrackMatchingPullPhi,
    kGMTrackMatchingPullTanl,
    kGMTrackMatchingPullInvQPt,
    kGMTracksP,
    kGMTrackDeltaTanl,
    kGMTrackDeltaPhi,
    kGMTrackDeltaPhiDeg,
    kGMTrackDeltaInvQPt,
    kGMTrackDeltaX,
    kGMTrackDeltaY,
    kGMTrackR,
    kGMTrackQ,
    kGMTrackRedChi2,
    kGMTrueTrackRedChi2,
    kGMFakeTrackRedChi2,
    kGMTrueTrackMatchChi2,
    kGMFakeTrackMatchChi2,
    kGMMatchChi2,
    kMCTrackspT,
    kMCTracksp,
    kMCTrackEta
  };

  std::map<int, const char*> TH1Names{
    {kGMTracksP, "GMTracks_P"},
    {kGMTrackPullX, "PullX"},
    {kGMTrackPullY, "PullY"},
    {kGMTrackPullPhi, "PullPhi"},
    {kGMTrackPullTanl, "PullTanl"},
    {kGMTrackPullInvQPt, "PullInvQPt"},
    {kGMTrackMatchingPullX, "MatchingPullX"},
    {kGMTrackMatchingPullY, "MatchingPullY"},
    {kGMTrackMatchingPullPhi, "MatchingPullPhi"},
    {kGMTrackMatchingPullTanl, "MatchingPullTanl"},
    {kGMTrackMatchingPullInvQPt, "MatchingPullInvQPt"},
    {kGMTrackDeltaTanl, "GMTracks_Delta_Tanl"},
    {kGMTrackDeltaPhi, "GMTracks_Delta_Phi"},
    {kGMTrackDeltaPhiDeg, "GMTracks_Delta_Phi_deg"},
    {kGMTrackDeltaInvQPt, "GMTracks_Delta_InvQPt"},
    {kGMTrackDeltaX, "GMTracks_Delta_X"},
    {kGMTrackDeltaY, "GMTracks_Delta_Y"},
    {kGMTrackR, "GMTracks_Delta_R"},
    {kGMTrackQ, "GMTracks_Charge_Match"},
    {kGMTrackRedChi2, "TracksReducedChi2"},
    {kGMTrueTrackRedChi2, "TracksReducedChi2_True"},
    {kGMFakeTrackRedChi2, "TracksReducedChi2_Fake"},
    {kGMTrueTrackMatchChi2, "MachingChi2_True"},
    {kGMFakeTrackMatchChi2, "MachingChi2_Fake"},
    {kGMMatchChi2, "MachingChi2"},
    {kMCTrackspT, "MCTracks_Pt"},
    {kMCTracksp, "MCTracks_P"},
    {kMCTrackEta, "MCTracks_eta"}};

  std::map<int, const char*> TH1Titles{
    {kGMTracksP, "Standalone Global Muon Tracks P"},
    {kGMTrackPullX, "PullX"},
    {kGMTrackPullY, "PullY"},
    {kGMTrackPullPhi, "PullPhi"},
    {kGMTrackPullTanl, "PullTanl"},
    {kGMTrackPullInvQPt, "PullInvQPt"},
    {kGMTrackMatchingPullX, "MatchingPullX"},
    {kGMTrackMatchingPullY, "MatchingPullY"},
    {kGMTrackMatchingPullPhi, "MatchingPullPhi"},
    {kGMTrackMatchingPullTanl, "MatchingPullTanl"},
    {kGMTrackMatchingPullInvQPt, "MatchingPullInvQPt"},
    {kGMTrackDeltaTanl, "tanl_{Fit} - tanl_{MC} "},
    {kGMTrackDeltaPhi, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaPhiDeg, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaInvQPt, "Global Muon Tracks \\Delta invQPt"},
    {kGMTrackDeltaX, "Global Muon Tracks Delta X at Z_vertex"},
    {kGMTrackDeltaY, "Global Muon Tracks Delta Y at Z_vertex"},
    {kGMTrackR, "Global Muon Tracks Delta R at Z_vertex"},
    {kGMTrackQ, "Global Muon Tracks Charge Match"},
    {kGMTrackRedChi2, "Global Muon Tracks ~ \\chi^2/DOF"},
    {kGMTrueTrackRedChi2, "Track \\chi^2/DOF ~ (true tracks) "},
    {kGMFakeTrackRedChi2, "Track \\chi^2/DOF ~ (fake tracks) "},
    {kGMTrueTrackMatchChi2, "Maching \\chi^2 ~ (true tracks) "},
    {kGMFakeTrackMatchChi2, "Maching \\chi^2 ~ (fake tracks) "},
    {kGMMatchChi2, "Global Muon Tracks MFT-MCH match ~ \\chi^2"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks Pseudorapidity"}};

  std::map<int, std::array<double, 3>> TH1Binning{
    {kGMTracksP, {500, pMin, pMax}},
    {kGMTrackPullX, {100, -5, 5}},
    {kGMTrackPullY, {100, -5, 5}},
    {kGMTrackPullPhi, {100, -5, +5}},
    {kGMTrackPullTanl, {100, -5, +5}},
    {kGMTrackPullInvQPt, {100, -5, +5}},
    {kGMTrackMatchingPullX, {100, -5, +5}},
    {kGMTrackMatchingPullY, {100, -5, +5}},
    {kGMTrackMatchingPullPhi, {100, -5, +5}},
    {kGMTrackMatchingPullTanl, {100, -5, +5}},
    {kGMTrackMatchingPullInvQPt, {100, -5, +5}},
    {kGMTrackDeltaTanl, {100, deltatanlMin, deltatanlMax}},
    {kGMTrackDeltaPhi, {100, deltaphiMin, deltaphiMax}},
    {kGMTrackDeltaPhiDeg, {100, TMath::RadToDeg() * deltaphiMin, TMath::RadToDeg() * deltaphiMax}},
    {kGMTrackDeltaInvQPt, {100, -1., 1.}},
    {kGMTrackDeltaX, {100, -.05, .05}},
    {kGMTrackDeltaY, {100, -.05, .05}},
    {kGMTrackR, {250, 0, 0.05}},
    {kGMTrackQ, {5, -2.1, 2.1}},
    {kGMTrackRedChi2, {1000, 0, 100}},
    {kGMTrueTrackRedChi2, {1000, 0, 100}},
    {kGMFakeTrackRedChi2, {1000, 0, 100}},
    {kGMTrueTrackMatchChi2, {1000, 0, 100}},
    {kGMFakeTrackMatchChi2, {1000, 0, 100}},
    {kGMMatchChi2, {1000, 0, 100}},
    {kMCTrackspT, {5000, 0, 50}},
    {kMCTracksp, {1000, pMin, pMax}},
    {kMCTrackEta, {1000, etaMin, etaMax}}};

  std::map<int, const char*> TH1XaxisTitles{
    {kGMTracksP, "p [GeV]"},
    {kGMTrackPullX, "Pull X"},
    {kGMTrackPullY, "Pull Y"},
    {kGMTrackPullPhi, "Pull Phi"},
    {kGMTrackPullTanl, "Pull Tanl"},
    {kGMTrackPullInvQPt, "Pull InvQPt"},
    {kGMTrackMatchingPullX, "Pull X"},
    {kGMTrackMatchingPullY, "Pull Y"},
    {kGMTrackMatchingPullPhi, "Pull Phi"},
    {kGMTrackMatchingPullTanl, "Pull Tanl"},
    {kGMTrackMatchingPullInvQPt, "Pull InvQPt"},
    {kGMTrackDeltaTanl, "\\Delta tanl"},
    {kGMTrackDeltaPhi, "\\Delta \\phi ~[rad]"},
    {kGMTrackDeltaPhiDeg, "\\Delta \\phi ~[deg]"},
    {kGMTrackDeltaInvQPt, "\\Delta invQPt"},
    {kGMTrackDeltaX, "\\Delta x ~[cm]"},
    {kGMTrackDeltaY, "\\Delta y ~[cm]"},
    {kGMTrackR, "\\Delta r ~[cm]"},
    {kGMTrackQ, "q_{fit}-q_{MC}"},
    {kGMTrackRedChi2, "\\chi^2/DOF"},
    {kGMTrueTrackRedChi2, "\\chi^2/DOF"},
    {kGMFakeTrackRedChi2, "\\chi^2/DOF"},
    {kGMTrueTrackMatchChi2, "\\chi^2"},
    {kGMFakeTrackMatchChi2, "\\chi^2"},
    {kGMMatchChi2, "\\chi^2"},
    {kMCTrackspT, "p_t [GeV]"},
    {kMCTracksp, "p [GeV]"},
    {kMCTrackEta, " \\eta"}};

  // Create histograms
  const int nTH1Histos = TH1Names.size();
  std::vector<std::unique_ptr<TH1F>> TH1Histos(nTH1Histos);
  auto nHisto = 0;
  for (auto& h : TH1Histos) {
    h = std::make_unique<TH1F>(TH1Names[nHisto], TH1Titles[nHisto],
                               (int)TH1Binning[nHisto][0],
                               TH1Binning[nHisto][1], TH1Binning[nHisto][2]);
    h->GetXaxis()->SetTitle(TH1XaxisTitles[nHisto]);
    ++nHisto;
  }

  const int nTH2Histos = TH2Names.size();
  std::vector<std::unique_ptr<TH2F>> TH2Histos(nTH2Histos);
  auto n2Histo = 0;
  for (auto& h : TH2Histos) {
    h = std::make_unique<TH2F>(TH2Names[n2Histo], TH2Titles[n2Histo],
                               (int)TH2Binning[n2Histo][0],
                               TH2Binning[n2Histo][1], TH2Binning[n2Histo][2],
                               (int)TH2Binning[n2Histo][3],
                               TH2Binning[n2Histo][4], TH2Binning[n2Histo][5]);
    // gStyle->SetLineWidth(4);
    // gROOT->ForceStyle();
    h->GetXaxis()->SetTitle(TH2XaxisTitles[n2Histo]);
    h->GetYaxis()->SetTitle(TH2YaxisTitles[n2Histo]);

    h->SetOption("COLZ");
    ++n2Histo;
  }

  // Profiles histograms
  auto PtRes_Profile = new TProfile("Pt_res_prof", "Profile of pt{fit}/pt{MC}",
                                    14, 0, 7, 0, 20, "s");
  PtRes_Profile->GetXaxis()->SetTitle("pt_{MC}");
  PtRes_Profile->GetYaxis()->SetTitle("mean(Pt_{Fit}/Pt_{MC})");

  auto DeltaX_Profile = new TProfile("DeltaX_prof", "Vertexing resolution", 14,
                                     0, 7, -1000., 1000., "s");
  DeltaX_Profile->GetXaxis()->SetTitle("pt_{MC} [GeV]");
  DeltaX_Profile->GetYaxis()->SetTitle("\\sigma_x ~[\\mu m]");

  // TEfficiency histogram
  TEfficiency* qMatchEff = new TEfficiency(
    "QMatchEff", "Charge Match;p_t [GeV];#epsilon", 40, 0, 20);

  TEfficiency* pairedMCHTracksEff = new TEfficiency(
    "PairingEff", "Paired_tracks;p_t [GeV];#epsilon", 40, 0, 20);
  TEfficiency* globalMuonCorrectMatchRatio =
    new TEfficiency("Correct_Match_Ratio",
                    " CorrectMatchRatio "
                    "(nCorrectMatches/NGlobalMuonTracks);p_t [GeV];#epsilon",
                    20, 0, 10);
  TEfficiency* globalMuonCombinedEff =
    new TEfficiency("Global_Matching_Efficiency",
                    "Global_Matching_Efficiency "
                    "(nCorrectMatches/NMCHTracks);p_t [GeV];#epsilon",
                    20, 0, 10);
  TEfficiency* closeMatchEff = new TEfficiency(
    "Close_Match_Eff", "Close Matches;p_t [GeV];#epsilon", 40, 0, 20);

  // Counters
  Int_t nChargeMatch = 0;
  Int_t nChargeMiss = 0;
  Int_t nCorrectMatchGMTracks = 0;
  Int_t nFakeGMTracks = 0;
  Int_t nNoMatchGMTracks = 0;

  // Files & Trees
  // MC
  TFile* o2sim_KineFileIn = new TFile(o2sim_KineFile.c_str());
  TTree* o2SimKineTree = (TTree*)o2sim_KineFileIn->Get("o2sim");
  TFile* sig_KineFileIn = new TFile(sig_KineFile.c_str());

  TTree* sigKineTree = (TTree*)sig_KineFileIn->Get("o2sim");

  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  vector<MCTrackT<float>>* mcTrSig = nullptr;
  sigKineTree->SetBranchAddress("MCTrack", &mcTrSig);

  Int_t numberOfEvents = o2SimKineTree->GetEntries();

  // Global Muon Tracks
  TFile* trkFileIn = new TFile(trkFile.c_str());
  TTree* gmTrackTree = (TTree*)trkFileIn->Get("GlobalFwdTracks");
  std::vector<GlobalMuonTrack> trackGMVec, *trackGMVecP = &trackGMVec;
  gmTrackTree->SetBranchAddress("fwdtracks", &trackGMVecP);

  vector<o2::MCCompLabel>* mcLabels = nullptr;
  gmTrackTree->SetBranchAddress("MCTruth", &mcLabels);

  std::string annotation = "";

  std::ifstream matcherConfig("MatchingConfig.txt");
  if (matcherConfig) {
    std::getline(matcherConfig, annotation);
    std::cout << "Matching config: " << annotation << std::endl;
  }
  matcherConfig.close();

  // MFT Tracks
  TFile* mfttrkFileIn = new TFile("mfttracks.root");
  TTree* mftTrackTree = (TTree*)mfttrkFileIn->Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackMFTVec, *trackMFTVecP = &trackMFTVec;
  mftTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);

  std::vector<o2::MCCompLabel>* mftMcLabels = nullptr;
  mftTrackTree->SetBranchAddress("MFTTrackMCTruth", &mftMcLabels);
  mftTrackTree->GetEntry(0);

  gmTrackTree->GetEntry(0);
  o2SimKineTree->GetEntry(0);
  sigKineTree->GetEntry(0);

  auto field_z = getZField(0, 0, -61.4); // Get field at Center of MFT

  std::string outfilename = "GlobalMuonChecks.root";
  TFile outFile(outfilename.c_str(), "RECREATE");

  // Reconstructed Global Muon Tracks
  std::cout << "Loop over Global Muon Tracks!" << std::endl;
  // GMTracks - Identify reconstructed tracks
  auto nCloseMatches = 0;
  auto iTrack = 0;

  if (0)
    for (auto& gmTrack : trackGMVec) {
      const auto& label = mcLabels->at(iTrack);
      std::cout << "iTrack = " << iTrack;
      label.print();
      iTrack++;
    }

  for (auto& gmTrack : trackGMVec) {

    auto bestMFTTrackMatchID = gmTrack.getMFTTrackID();
    auto& mftTrackMatch = trackMFTVec[bestMFTTrackMatchID];
    MCTrackT<float>* thisTrack;
    auto trackRedChi2 = gmTrack.getTrackChi2() / (2 * mftTrackMatch.getNumberOfPoints() - 5);
    auto matchChi2 = gmTrack.getMatchingChi2();
    TH1Histos[kGMTrackRedChi2]->Fill(trackRedChi2);
    TH1Histos[kGMMatchChi2]->Fill(matchChi2);

    const auto& label = mcLabels->at(iTrack);

    // std::cout << "iTrack = " << iTrack;
    // label.print();

    //if (iEvent == label.getEventID()) {
    if (DEBUG_VERBOSE) {
      // std::cout << "  Global Track ID = " <<  iTrack << " ; MFTMatchID =
      // " << bestMFTTrackMatchID << " SourceID = " <<
      // label.getSourceID()
      // << " ; EventID = " << label.getEventID() << ":  trackID = " <<
      // label.getTrackID() << " ; isFake = " << label.isFake() << "
      // Label: ";
      std::cout << "  Global Track ID = " << iTrack
                << " ; MFTMatchID = " << bestMFTTrackMatchID << " Label: ";
      label.print();

      // std::cout << "        bestMFTTrackMatchID = " <<
      // bestMFTTrackMatchID << " / labelMFTBestMatch = ";
      // labelMFTBestMatch.print();
    }
    if (gmTrack.isCloseMatch())
      nCloseMatches++;

    auto thisTrkID = label.getTrackID();
    auto thisEvtID = label.getEventID();
    if (label.getSourceID() == 0) {
      o2SimKineTree->GetEntry(thisEvtID);
      thisTrack = &(mcTr->at(thisTrkID));
    } else if (label.getSourceID() == 1) {
      sigKineTree->GetEntry(thisEvtID);
      thisTrack = &(mcTrSig->at(thisTrkID));
    } else {
      iTrack++;
      continue;
    }

    pairedMCHTracksEff->Fill(bestMFTTrackMatchID > -1, gmTrack.getPt());
    globalMuonCombinedEff->Fill(label.isCorrect(), gmTrack.getPt());
    closeMatchEff->Fill(gmTrack.isCloseMatch(), gmTrack.getPt());

    if (bestMFTTrackMatchID >= 0) {
      globalMuonCorrectMatchRatio->Fill(label.isCorrect(),
                                        gmTrack.getPt());
    }
    if (label.isCorrect()) { // Correct match track: add to histograms
      nCorrectMatchGMTracks++;

      auto vx_MC = thisTrack->GetStartVertexCoordinatesX();
      auto vy_MC = thisTrack->GetStartVertexCoordinatesY();
      auto vz_MC = thisTrack->GetStartVertexCoordinatesZ();
      auto Pt_MC = thisTrack->GetPt();
      auto P_MC = thisTrack->GetP();
      auto phi_MC = TMath::ATan2(thisTrack->Py(), thisTrack->Px());
      auto eta_MC =
        atanh(thisTrack->GetStartVertexMomentumZ() / P_MC); // eta;
      auto tanl_MC = thisTrack->Pz() / thisTrack->GetPt();
      auto pdgcode_MC = thisTrack->GetPdgCode();
      // std::cout << "pdgcode_MC = " <<  pdgcode_MC;
      int Q_MC;
      if (TDatabasePDG::Instance()->GetParticle(pdgcode_MC)) {
        Q_MC =
          TDatabasePDG::Instance()->GetParticle(pdgcode_MC)->Charge() / 3;
        if (DEBUG_VERBOSE)
          std::cout << "      => "
                    << TDatabasePDG::Instance()
                         ->GetParticle(pdgcode_MC)
                         ->GetName()
                    << "\n";
      } else {
        Q_MC = 0;
        std::cout << " => pdgcode ERROR " << Q_MC << "\n";
      }

      gmTrack.propagateToZhelix(vz_MC, field_z);
      // gmTrack.propagateToZquadratic(vz_MC,field_z);
      // gmTrack.propagateToZlinear(vz_MC,field_z);

      auto Q_fit = gmTrack.getCharge();
      auto dx = gmTrack.getX() - vx_MC;
      auto dy = gmTrack.getY() - vy_MC;
      auto d_eta = gmTrack.getEta() - eta_MC;
      auto d_tanl = gmTrack.getTanl() - tanl_MC;
      auto Pt_fit = gmTrack.getPt();
      auto d_invQPt = Q_fit / Pt_fit - Q_MC / Pt_MC;
      auto P_fit = gmTrack.getP();
      auto P_res = P_fit / P_MC;
      auto Pt_res = Pt_fit / Pt_MC;
      auto d_Phi = gmTrack.getPhi() - phi_MC;
      auto d_Charge = Q_fit - Q_MC;
      TH1Histos[kGMTracksP]->Fill(gmTrack.getP());
      TH1Histos[kGMTrackDeltaTanl]->Fill(d_tanl);
      TH1Histos[kGMTrackDeltaPhi]->Fill(d_Phi);
      TH1Histos[kGMTrackDeltaInvQPt]->Fill(d_invQPt);
      TH1Histos[kGMTrackDeltaPhiDeg]->Fill(TMath::RadToDeg() * d_Phi);
      TH1Histos[kGMTrackDeltaX]->Fill(dx);

      TH1Histos[kGMTrackPullX]->Fill(dx / sqrt(gmTrack.getCovariances()(0, 0)));
      TH1Histos[kGMTrackPullY]->Fill(dy / sqrt(gmTrack.getCovariances()(1, 1)));
      TH1Histos[kGMTrackPullPhi]->Fill(d_Phi / sqrt(gmTrack.getCovariances()(2, 2)));
      TH1Histos[kGMTrackPullTanl]->Fill(d_tanl / sqrt(gmTrack.getCovariances()(3, 3)));
      TH1Histos[kGMTrackPullInvQPt]->Fill(d_invQPt / sqrt(gmTrack.getCovariances()(4, 4)));

      //
      auto Residuals2Cov = gmTrack.computeResiduals2Cov(trackMFTVec[bestMFTTrackMatchID].getOutParam());
      TH1Histos[kGMTrackMatchingPullX]->Fill(Residuals2Cov(0));
      TH1Histos[kGMTrackMatchingPullY]->Fill(Residuals2Cov(1));
      TH1Histos[kGMTrackMatchingPullPhi]->Fill(Residuals2Cov(2));
      TH1Histos[kGMTrackMatchingPullTanl]->Fill(Residuals2Cov(3));
      TH1Histos[kGMTrackMatchingPullInvQPt]->Fill(Residuals2Cov(4));

      DeltaX_Profile->Fill(Pt_MC, dx * 1e4);
      TH1Histos[kGMTrackDeltaY]->Fill(dy);
      TH1Histos[kGMTrackR]->Fill(sqrt(dx * dx + dy * dy));
      TH1Histos[kGMTrackQ]->Fill(d_Charge);
      TH1Histos[kGMTrueTrackRedChi2]->Fill(trackRedChi2);
      TH1Histos[kGMTrueTrackMatchChi2]->Fill(matchChi2);
      TH2Histos[kGMTrackDeltaXYVertex]->Fill(10. * dx, 10. * dy);
      TH2Histos[kGMTrackQPtRec_MC]->Fill(Pt_MC * Q_MC, Pt_fit * Q_fit);
      TH2Histos[kGMTrackPtResolution]->Fill(Pt_MC, Pt_fit / Pt_MC);
      PtRes_Profile->Fill(Pt_MC, Pt_fit / Pt_MC);
      TH2Histos[kGMTrackInvPtResolution]->Fill(
        Pt_MC, (1.0 / Pt_fit - 1.0 / Pt_MC) * Pt_MC);

      // MC histos
      TH1Histos[kMCTrackspT]->Fill(Pt_MC);
      TH1Histos[kMCTracksp]->Fill(P_MC);
      TH1Histos[kMCTrackEta]->Fill(eta_MC);
      TH2Histos[kMCTracksEtaZ]->Fill(vz_MC, eta_MC);

      d_Charge ? nChargeMiss++ : nChargeMatch++;
      qMatchEff->Fill(!d_Charge, Pt_MC);
    } else {
      if (bestMFTTrackMatchID >= 0) {
        nFakeGMTracks++;
        TH1Histos[kGMFakeTrackMatchChi2]->Fill(matchChi2);
        TH1Histos[kGMFakeTrackRedChi2]->Fill(trackRedChi2);
      } else
        nNoMatchGMTracks++;
    }
    //}
    iTrack++;

  } // Loop on GMTracks
  // }   // Loop over events

  Int_t nRecoGMTracks = nCorrectMatchGMTracks + nFakeGMTracks;
  Int_t nMCHTracks = nRecoGMTracks + nNoMatchGMTracks;

  // Customize histograms
  TH1Histos[kGMTrackQ]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch,
         100. * nChargeMatch / (nChargeMiss + nChargeMatch)));

  qMatchEff->SetTitle(Form("Charge match = %.2f%%",
                           100. * nChargeMatch / (nChargeMiss + nChargeMatch)));
  pairedMCHTracksEff->SetTitle(
    Form("Paired_MCH_tracks_=_%.2f%%", 100. * nRecoGMTracks / (nMCHTracks)));
  globalMuonCorrectMatchRatio->SetTitle(
    Form("Correct_Match_Ratio = %.2f%%",
         100. * nCorrectMatchGMTracks / (nRecoGMTracks)));
  closeMatchEff->SetTitle(
    Form("Close_Match_=_%.2f%%", 100. * nCloseMatches / (nMCHTracks)));

  // Remove stat boxes
  TH2Histos[kGMTrackQPtRec_MC]->SetStats(0);
  TH2Histos[kGMTrackPtResolution]->SetStats(0);
  TH2Histos[kGMTrackInvPtResolution]->SetStats(0);
  TH2Histos[kMCTracksEtaZ]->SetStats(0);
  PtRes_Profile->SetStats(0);
  DeltaX_Profile->SetStats(0);
  TH1Histos[kGMTrackQ]->SetStats(0);

  // Fit Slices: Pt resolution
  FitSlicesy(*TH2Histos[kGMTrackInvPtResolution], *TH2Histos[kGMTrackQPtRec_MC]);
  FitSlicesy(*TH2Histos[kGMTrackPtResolution], *TH2Histos[kGMTrackQPtRec_MC]);

  // sigmaX resultion Profile
  TH1D* DeltaX_Error;
  DeltaX_Error = DeltaX_Profile->ProjectionX("DeltaX_Error", "C=E");
  DeltaX_Error->GetYaxis()->SetTitleOffset(1.25);
  DeltaX_Error->SetMaximum(500);
  // DeltaX_Error->GetYaxis()->SetLimits(0,500.0);

  // Summary Canvases
  // Matching summary
  auto matching_summary = summary_report_3x2(
    *pairedMCHTracksEff, *globalMuonCorrectMatchRatio, *closeMatchEff,
    *TH2Histos[kGMTrackDeltaXYVertex], *DeltaX_Error, *PtRes_Profile,
    "MatchingSummary", annotation, 0, 0, 0, 0, 0, 0, "-", "-", "-",
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex]->GetEntries()),
    "-", "-");

  // Parameters resolution
  auto param_resolution = summary_report_3x2(
    *TH2Histos[kGMTrackDeltaXYVertex], *TH2Histos[kGMTrackPtResolution],
    *PtRes_Profile, *DeltaX_Error, *TH2Histos[kGMTrackQPtRec_MC], *qMatchEff,
    "ParamSummary", annotation, 0, 0, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackPtResolution]->Integral() /
                     TH2Histos[kGMTrackPtResolution]->GetEntries()),
    "-", "-",
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPtRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPtRec_MC]->GetEntries()),
    "-");

  // Covariances summary
  auto covariances_summary = summary_report_3x2(
    *TH1Histos[kGMTrackPullX], *TH1Histos[kGMTrackPullPhi],
    *TH1Histos[kGMTrackPullInvQPt], *TH1Histos[kGMTrackPullY],
    *TH1Histos[kGMTrackPullTanl], *TH2Histos[kGMTrackQPtRec_MC],
    "CovariancesSummary", annotation, 1, 1, 1, 1, 1, 0,
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullX]->Integral() /
                     TH1Histos[kGMTrackPullX]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullPhi]->Integral() /
                     TH1Histos[kGMTrackPullPhi]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullInvQPt]->Integral() /
                     TH1Histos[kGMTrackPullInvQPt]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullY]->Integral() /
                     TH1Histos[kGMTrackPullY]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullTanl]->Integral() /
                     TH1Histos[kGMTrackPullTanl]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPtRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPtRec_MC]->GetEntries()));

  // MCH Residuals Covariances summary
  auto MCHcovariances_summary = summary_report_3x2(
    *TH1Histos[kGMTrackMatchingPullX], *TH1Histos[kGMTrackMatchingPullPhi],
    *TH1Histos[kGMTrackMatchingPullInvQPt], *TH1Histos[kGMTrackMatchingPullY],
    *TH1Histos[kGMTrackMatchingPullTanl], *TH2Histos[kGMTrackQPtRec_MC],
    "MatchingPullsSummary", annotation, 1, 1, 1, 1, 1, 0,
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackMatchingPullX]->Integral() /
                     TH1Histos[kGMTrackMatchingPullX]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackMatchingPullPhi]->Integral() /
                     TH1Histos[kGMTrackMatchingPullPhi]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackMatchingPullInvQPt]->Integral() /
                     TH1Histos[kGMTrackMatchingPullInvQPt]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackMatchingPullY]->Integral() /
                     TH1Histos[kGMTrackMatchingPullY]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackMatchingPullTanl]->Integral() /
                     TH1Histos[kGMTrackMatchingPullTanl]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPtRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPtRec_MC]->GetEntries()));

  // Covariances summary 3x3
  auto par_cov_summary3x3 = summary_report_3x3(
    *TH2Histos[kGMTrackDeltaXYVertex], *TH1Histos[kGMTrackPullX],
    *TH1Histos[kGMTrackPullY], *DeltaX_Error,
    *TH2Histos[kGMTrackQPtRec_MC], *TH1Histos[kGMTrackPullPhi], *qMatchEff,
    *TH1Histos[kGMTrackPullInvQPt], *TH1Histos[kGMTrackPullTanl],
    "ParametersAndCovariancesSummary_3x3", annotation, 0, 1, 1, 0, 0, 1, 0, 1, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullX]->Integral() /
                     TH1Histos[kGMTrackPullX]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullY]->Integral() /
                     TH1Histos[kGMTrackPullY]->GetEntries()),
    "-",
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPtRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPtRec_MC]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullPhi]->Integral() /
                     TH1Histos[kGMTrackPullPhi]->GetEntries()),
    "-",
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullInvQPt]->Integral() /
                     TH1Histos[kGMTrackPullInvQPt]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPullTanl]->Integral() /
                     TH1Histos[kGMTrackPullTanl]->GetEntries()));

  auto pt_resolution = summary_report(
    *TH2Histos[kGMTrackPtResolution], *TH2Histos[kGMTrackQPtRec_MC],
    *PtRes_Profile, *qMatchEff, "PtSummary", annotation, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackPtResolution]->Integral() /
                     TH2Histos[kGMTrackPtResolution]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPtRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPtRec_MC]->GetEntries()));

  auto invpt_resolution = summary_report(
    *TH2Histos[kGMTrackInvPtResolution], *TH2Histos[kGMTrackQPtRec_MC],
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kGMTrackInvPtResolution]->GetName()) +
       std::string("_1"))
        .c_str()),
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kGMTrackInvPtResolution]->GetName()) +
       std::string("_2"))
        .c_str()),
    "InvPtSummary", annotation, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackInvPtResolution]->Integral() /
                     TH2Histos[kGMTrackInvPtResolution]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPtRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPtRec_MC]->GetEntries()));

  auto vertexing_resolution = summary_report(
    *TH2Histos[kGMTrackDeltaXYVertex], *TH1Histos[kGMTrackDeltaX],
    *DeltaX_Error, *TH1Histos[kGMTrackDeltaPhiDeg], "VertexingSummary",
    annotation, 0, 1, 0, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaX]->Integral() /
                     TH1Histos[kGMTrackDeltaX]->GetEntries()),
    Form("-"),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiDeg]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiDeg]->GetEntries()));

  // Write histograms to file and export images

  outFile.mkdir("MoreHistos");
  outFile.cd("MoreHistos");

  for (auto& h : TH2Histos) {
    h->Write();
    if (EXPORT_HISTOS_IMAGES)
      exportHisto(*h);
  }

  for (auto& h : TH1Histos) {
    h->Write();
    if (EXPORT_HISTOS_IMAGES)
      exportHisto(*h);
  }

  PtRes_Profile->Write();
  DeltaX_Profile->Write();
  DeltaX_Error->Write();
  qMatchEff->Write();
  pairedMCHTracksEff->Write();
  globalMuonCorrectMatchRatio->Write();
  closeMatchEff->Write();
  globalMuonCombinedEff->Write();
  outFile.cd();
  //outFile.WriteObjectAny(&matching_helper, "MatchingHelper", "Matching Helper");

  outFile.Close();

  std::cout << std::endl;
  std::cout << "---------------------------------------------------"
            << std::endl;
  std::cout << "-------------   Matching Summary   ----------------"
            << std::endl;
  std::cout << "---------------------------------------------------"
            << std::endl;
  std::cout << " P_mean = " << TH1Histos[kGMTracksP]->GetMean() << std::endl;
  std::cout << " P_StdDev = " << TH1Histos[kGMTracksP]->GetStdDev()
            << std::endl;
  std::cout << " Tanl_mean = " << TH1Histos[kGMTrackDeltaTanl]->GetMean()
            << std::endl;
  std::cout << " Tanl_StdDev = " << TH1Histos[kGMTrackDeltaTanl]->GetStdDev()
            << std::endl;
  std::cout << " Phi_mean = " << TH1Histos[kGMTrackDeltaPhi]->GetMean()
            << std::endl;
  std::cout << " Phi_StdDev = " << TH1Histos[kGMTrackDeltaPhi]->GetStdDev()
            << std::endl;
  std::cout << " Phi_meanDeg = " << TH1Histos[kGMTrackDeltaPhiDeg]->GetMean()
            << std::endl;
  std::cout << " Phi_StdDevDeg = "
            << TH1Histos[kGMTrackDeltaPhiDeg]->GetStdDev() << std::endl;
  std::cout << " DeltaX_mean = " << TH1Histos[kGMTrackDeltaX]->GetMean()
            << std::endl;
  std::cout << " DeltaX_StdDev = " << TH1Histos[kGMTrackDeltaX]->GetStdDev()
            << std::endl;
  std::cout << " DeltaY_mean = " << TH1Histos[kGMTrackDeltaY]->GetMean()
            << std::endl;
  std::cout << " DeltaY_StdDev = " << TH1Histos[kGMTrackDeltaY]->GetStdDev()
            << std::endl;
  std::cout << " R_mean = " << TH1Histos[kGMTrackR]->GetMean() << std::endl;
  std::cout << " R_StdDev = " << TH1Histos[kGMTrackR]->GetStdDev() << std::endl;
  std::cout << " Charge_mean = " << TH1Histos[kGMTrackDeltaY]->GetMean()
            << std::endl;
  std::cout << " nChargeMatch = " << nChargeMatch << " ("
            << 100. * nChargeMatch / (nChargeMiss + nChargeMatch) << "%)"
            << std::endl;
  std::cout << " nTrackMatch = " << nCorrectMatchGMTracks << " ("
            << 100. * nCorrectMatchGMTracks / (nRecoGMTracks) << "%)"
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "------------"
            << std::endl;

  std::cout << std::endl;
  std::cout << "---------------------------------------------------------------"
               "------------"
            << std::endl;
  std::cout << "------------------------   Track matching Summary   "
               "-----------------------"
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "------------"
            << std::endl;
  std::cout << " ==> " << nMCHTracks << " MCH Tracks in " << numberOfEvents
            << " events" << std::endl;
  std::cout << " ==> " << nNoMatchGMTracks
            << " dangling MCH Tracks (no MFT track to match)"
            << " (" << 100. * nNoMatchGMTracks / (nMCHTracks) << "%)"
            << std::endl;
  std::cout << " ==> " << nRecoGMTracks << " reconstructed Global Muon Tracks"
            << " (" << 100. * nRecoGMTracks / (nMCHTracks) << "%)" << std::endl;
  std::cout << " ==> " << nFakeGMTracks << " fake Global Muon Tracks"
            << " (contamination = " << 100. * nFakeGMTracks / (nRecoGMTracks)
            << "%)" << std::endl;
  std::cout << " ==> " << nCloseMatches
            << " close matches - correct MFT track in search window"
            << " (" << 100. * nCloseMatches / (nMCHTracks) << "%)"
            << std::endl;
  std::cout << " ==> " << nCorrectMatchGMTracks
            << " Correct Match Global Muon Tracks"
            << " (Correct_Match_Ratio = "
            << 100. * nCorrectMatchGMTracks / (nRecoGMTracks) << "%)"
            << " (eff. = " << 100. * nCorrectMatchGMTracks / (nMCHTracks)
            << "%)" << std::endl;

  std::cout << "---------------------------------------------------------------"
               "-----------"
            << std::endl;
  std::cout << " Annotation: " << annotation << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------"
            << std::endl;
  std::cout << std::endl;

  /*
  std::cout << "matching_helper.nMCHTracks = " << matching_helper.nMCHTracks <<
  std::endl; std::cout << "matching_helper.nNoMatch = " <<
  matching_helper.nNoMatch << std::endl; std::cout <<
  "matching_helper.nGMTracks() = " << matching_helper.nGMTracks() << std::endl;
  std::cout << "matching_helper.nFakes = " << matching_helper.nFakes <<
  std::endl; std::cout << "matching_helper.nCorrectMatches = " <<
  matching_helper.nCorrectMatches << std::endl; std::cout <<
  "matching_helper.nCloseMatches = " << matching_helper.nCloseMatches <<
  std::endl; std::cout << "matching_helper.getCorrectMatchRatio() = " <<
  matching_helper.getCorrectMatchRatio() << std::endl; std::cout <<
  "matching_helper.getPairingEfficiency() = " <<
  matching_helper.getPairingEfficiency() << std::endl;
  */

  return 0;
}
