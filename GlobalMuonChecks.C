#if !defined(__CLING__) || defined(__ROOTCLING__)

#ifdef __MAKECINT__
#pragma link C++ class GlobalMuonTrack + ;
#pragma link C++ class std::vector < GlobalMuonTrack> + ;
#pragma link C++ class MatchingHelper + ;
#endif

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
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

#include "include/GlobalMuonTrack.h"
#include "macrohelpers/HistosHelpers.C"
#include "macrohelpers/MagField.C"

using o2::MCTrackT;
using GlobalMuonTrack = o2::track::GlobalMuonTrack;
using eventFoundTracks = std::vector<bool>;
using std::vector;
vector<eventFoundTracks> allFoundGMTracks; // True for reconstructed tracks -
                                           // one vector of bool per event

bool DEBUG_VERBOSE = false;
bool EXPORT_HISTOS_IMAGES = false;

//_________________________________________________________________________________________________
int GlobalMuonChecks(const std::string trkFile = "GlobalMuonTracks.root",
                     const std::string o2sim_KineFile = "o2sim_Kine.root")
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
  Double_t deltaphiMin = -.2; //-3.15,
  Double_t deltaphiMax = .2;  //+3.15,
  Double_t deltatanlMin = -2.0;
  Double_t deltatanlMax = 2.0;

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

  enum TH2HistosCodes {
    kGMTrackDeltaXYVertex,
    kGMTrackDeltaXYVertex0_1,
    kGMTrackDeltaXYVertex1_4,
    kGMTrackDeltaXYVertex4plus,
    kGMTrackChi2vsFitChi2,
    kGMTrackQPRec_MC,
    kGMTrackPtResolution,
    kGMTrackInvPtResolution,
    kMCTracksEtaZ
  };

  std::map<int, const char*> TH2Names{
    {kGMTrackDeltaXYVertex, "Global Muon Tracks Vertex at Z = 0"},
    {kGMTrackDeltaXYVertex0_1, "Global Muon Tracks Vertex at Z = 0 Pt0_1"},
    {kGMTrackDeltaXYVertex1_4, "Global Muon Tracks Vertex at Z = 0 Pt1_4"},
    {kGMTrackDeltaXYVertex4plus,
     "Global Muon Tracks Vertex at Z = 0 Pt4plus"},
    {kGMTrackChi2vsFitChi2, "Global Muon TracksChi2vsFitChi2"},
    {kGMTrackQPRec_MC, "GM Track QP FITxMC"},
    {kGMTrackPtResolution, "GM Track Pt Resolution"},
    {kGMTrackInvPtResolution, "GM Track InvPt Resolution"},
    {kMCTracksEtaZ, "MCTracks_eta_z"}};

  std::map<int, const char*> TH2Titles{
    {kGMTrackDeltaXYVertex, "Global Muon Tracks at Z_vertex"},
    {kGMTrackDeltaXYVertex0_1, "Global Muon Tracks at Z_vertex (pt < 1)"},
    {kGMTrackDeltaXYVertex1_4, "Global Muon Tracks at Z_vertex (1 < pt < 4)"},
    {kGMTrackDeltaXYVertex4plus, "Global Muon Tracks at Z_vertex (pt > 4)"},
    {kGMTrackChi2vsFitChi2, "Tracks Chi2 vs FitChi2"},
    {kGMTrackQPRec_MC, "Charged Momentum: Reconstructed vs MC"},
    {kGMTrackPtResolution, "Pt Resolution"},
    {kGMTrackInvPtResolution, "InvPt Resolution"},
    {kMCTracksEtaZ, "MC Tracks: Pseudorapidity vs zVertex"}};

  std::map<int, std::array<double, 6>> TH2Binning{
    {kGMTrackDeltaXYVertex, {100, -.5, .5, 100, -.5, .5}},
    {kGMTrackDeltaXYVertex0_1, {100, -.5, .5, 100, -.5, .5}},
    {kGMTrackDeltaXYVertex1_4, {100, -.5, .5, 100, -.5, .5}},
    {kGMTrackDeltaXYVertex4plus, {100, -.5, .5, 100, -.5, .5}},
    {kGMTrackChi2vsFitChi2, {500, 0, 1000, 250, 0., 500.}},
    {kGMTrackQPRec_MC, {50, -100, 100, 50, -100, 100}},
    {kGMTrackPtResolution, {20, 0, 10, 100, 0, 5}},
    {kGMTrackInvPtResolution, {14, 0, 7, 300, -2, 2}},
    {kMCTracksEtaZ, {31, -15, 16, 25, etaMin, etaMax}}};

  std::map<int, const char*> TH2XaxisTitles{
    {kGMTrackDeltaXYVertex, "\\Delta x ~[mm]"},
    {kGMTrackDeltaXYVertex0_1, "\\Delta x ~[mm]"},
    {kGMTrackDeltaXYVertex1_4, "\\Delta x ~[mm]"},
    {kGMTrackDeltaXYVertex4plus, "\\Delta x ~[mm]"},
    {kGMTrackChi2vsFitChi2, "Fit ~ \\chi^2"},
    {kGMTrackQPRec_MC, "(q.p)_{MC} [GeV]"},
    {kGMTrackPtResolution, "pt_{MC} [GeV]"},
    {kGMTrackInvPtResolution, "pt_{MC} [GeV]"},
    {kMCTracksEtaZ, "Vertex PosZ [cm]"}};

  std::map<int, const char*> TH2YaxisTitles{
    {kGMTrackDeltaXYVertex, "\\Delta y ~[mm]"},
    {kGMTrackDeltaXYVertex0_1, "\\Delta y ~[mm]"},
    {kGMTrackDeltaXYVertex1_4, "\\Delta y ~[mm]"},
    {kGMTrackDeltaXYVertex4plus, "\\Delta y ~[mm]"},
    {kGMTrackChi2vsFitChi2, "Track ~ \\chi^2"},
    {kGMTrackQPRec_MC, "(q.p)_{fit} [GeV]"},
    {kGMTrackPtResolution, "pt_{fit} / pt_{MC}"},
    {kGMTrackInvPtResolution, "(1/(p_t)_{fit} - 1/(p_t)_{MC})*(p_t)_{MC}"},
    {kMCTracksEtaZ, "\\eta"}};

  enum TH1HistosCodes {
    kGMTrackDeltaXErr,
    kGMTrackDeltaYErr,
    kGMTrackDeltaPhiErr,
    kGMTrackDeltaTanLErr,
    kGMTrackDeltainvQPtErr,
    kMCHResTrackDeltaXErr,
    kMCHResTrackDeltaYErr,
    kMCHResTrackDeltaPhiErr,
    kMCHResTrackDeltaTanLErr,
    kMCHResTrackDeltainvQPtErr,
    kGMTrackXChi2,
    kGMTrackYChi2,
    kGMTrackPhiChi2,
    kGMTrackTanlChi2,
    kGMTrackinvQPtChi2,
    kFitChi2,
    kGMTracksP,
    kGMTrackDeltaTanl,
    kGMTrackDeltaTanl0_1,
    kGMTrackDeltaTanl1_4,
    kGMTrackDeltaTanl4plus,
    kGMTrackDeltaPhi,
    kGMTrackDeltaPhi0_1,
    kGMTrackDeltaPhi1_4,
    kGMTrackDeltaPhi4plus,
    kGMTrackDeltaPhiDeg,
    kGMTrackDeltaPhiDeg0_1,
    kGMTrackDeltaPhiDeg1_4,
    kGMTrackDeltaPhiDeg4plus,
    kGMTrackDeltaInvQPt,
    kGMTrackDeltaX,
    kGMTrackDeltaX0_1,
    kGMTrackDeltaX1_4,
    kGMTrackDeltaX4plus,
    kGMTrackDeltaY,
    kGMTrackR,
    kGMTrackQ,
    kGMTrackQ0_1,
    kGMTrackQ1_4,
    kGMTrackQ4plus,
    kGMTrackChi2,
    kMCTrackspT,
    kMCTracksp,
    kMCTrackEta
  };

  std::map<int, const char*> TH1Names{
    {kGMTracksP, "Global Muon Tracks Fitted p"},
    {kGMTrackDeltaXErr, "Delta X / SigmaX"},
    {kGMTrackDeltaYErr, "Delta Y / SigmaY"},
    {kGMTrackDeltaPhiErr, "Delta Phi at Vertex / SigmaPhi"},
    {kGMTrackDeltaTanLErr, "Delta_Tanl / SigmaTanl"},
    {kGMTrackDeltainvQPtErr, "Delta_InvQPt / Sigma_{q/pt}"},
    {kMCHResTrackDeltaXErr, "MCH Delta X / SigmaX"},
    {kMCHResTrackDeltaYErr, "MCH Delta Y / SigmaY"},
    {kMCHResTrackDeltaPhiErr, "MCH Delta Phi at Vertex / SigmaPhi"},
    {kMCHResTrackDeltaTanLErr, "MCH Delta_Tanl / SigmaTanl"},
    {kMCHResTrackDeltainvQPtErr, "MCH Delta_InvQPt / Sigma_{q/pt}"},
    {kGMTrackDeltaTanl, "Global Muon Tracks Fitted Delta_tanl"},
    {kGMTrackXChi2, "X Chi2"},
    {kGMTrackYChi2, "Y Chi2"},
    {kGMTrackPhiChi2, "Phi chi2"},
    {kGMTrackTanlChi2, "Tanl Chi2"},
    {kGMTrackinvQPtChi2, "InvQPt Chi2"},
    {kFitChi2, "Fit Chi2"},
    {kGMTrackDeltaTanl0_1, "Global Muon Tracks tanl (pt < 1)"},
    {kGMTrackDeltaTanl1_4, "Global Muon Tracks tanl (1 < pt < 4)"},
    {kGMTrackDeltaTanl4plus, "Global Muon Tracks tanl (pt > 4)"},
    {kGMTrackDeltaPhi, "Global Muon Tracks Fitted Phi at Vertex"},
    {kGMTrackDeltaPhi0_1,
     "Global Muon Tracks Fitted Phi at Vertex [rad] (pt < 1)"},
    {kGMTrackDeltaPhi1_4,
     "Global Muon Tracks Fitted Phi at Vertex [rad] (1 < pt < 4)"},
    {kGMTrackDeltaPhi4plus,
     "Global Muon Tracks Fitted Phi at Vertex [rad] (pt > 4)"},
    {kGMTrackDeltaPhiDeg, "Global Muon Tracks Fitted Phi at Vertex [deg]"},
    {kGMTrackDeltaPhiDeg0_1,
     "Global Muon Tracks Fitted Phi at Vertex [deg] (pt < 1)"},
    {kGMTrackDeltaPhiDeg1_4,
     "Global Muon Tracks Fitted Phi at Vertex [deg] (1 < pt < 4)"},
    {kGMTrackDeltaPhiDeg4plus,
     "Global Muon Tracks Fitted Phi at Vertex [deg] (pt > 4)"},
    {kGMTrackDeltaInvQPt, "Global Muon Tracks invQPt"},
    {kGMTrackDeltaX, "Global Muon Tracks Delta X"},
    {kGMTrackDeltaX0_1, "Global Muon Tracks Delta X (pt < 1)"},
    {kGMTrackDeltaX1_4, "Global Muon Tracks Delta X (1 < pt < 4)"},
    {kGMTrackDeltaX4plus, "Global Muon Tracks Delta X (pt > 4)"},
    {kGMTrackDeltaY, "Global Muon Tracks Delta Y"},
    {kGMTrackR, "Global Muon Tracks Delta R"},
    {kGMTrackQ, "Charge Match"},
    {kGMTrackQ0_1, "Charge Match (pt < 1)"},
    {kGMTrackQ1_4, "Charge Match (1 < pt < 4)"},
    {kGMTrackQ4plus, "Charge Match (pt > 4)"},
    {kGMTrackChi2, "Tracks Chi2"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks eta"}};

  std::map<int, const char*> TH1Titles{
    {kGMTracksP, "Standalone Global Muon Tracks P"},
    {kGMTrackDeltaXErr, "\\Delta X / \\sigma_X"},
    {kGMTrackDeltaYErr, "\\Delta Y / \\sigma_Y"},
    {kGMTrackDeltaPhiErr, "\\Delta \\phi / \\sigma_\\phi"},
    {kGMTrackDeltaTanLErr, "\\Delta TanL / \\sigma_{TanL} "},
    {kGMTrackDeltainvQPtErr, "\\Delta(q/Pt) / \\sigma_{q/pt}"},
    {kMCHResTrackDeltaXErr, "\\Delta X / \\sigma_X"},
    {kMCHResTrackDeltaYErr, "\\Delta Y / \\sigma_Y"},
    {kMCHResTrackDeltaPhiErr, "\\Delta \\phi / \\sigma_\\phi"},
    {kMCHResTrackDeltaTanLErr, "\\Delta TanL / \\sigma_{TanL} "},
    {kMCHResTrackDeltainvQPtErr, "\\Delta(q/Pt) / \\sigma_{q/pt}"},
    {kGMTrackXChi2, "\\chi^2(x)"},
    {kGMTrackYChi2, "\\chi^2(y)"},
    {kGMTrackPhiChi2, "\\chi^2(\\phi)"},
    {kGMTrackTanlChi2, "\\chi^2(TanL)"},
    {kGMTrackinvQPtChi2, "\\chi^2(InvQP_t)"},
    {kFitChi2, "Fit Chi2"},
    {kGMTrackDeltaTanl, "tanl_{Fit} - tanl_{MC} "},
    {kGMTrackDeltaTanl0_1, "tanl_{Fit} - tanl_{MC} (pt < 1)"},
    {kGMTrackDeltaTanl1_4, "tanl_{Fit} - tanl_{MC} (1 < p_t < 4)"},
    {kGMTrackDeltaTanl4plus, "tanl_{Fit} - tanl_{MC} (p_t > 4)"},
    {kGMTrackDeltaPhi, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaPhi0_1, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaPhi1_4, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaPhi4plus, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaPhiDeg, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaPhiDeg0_1, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaPhiDeg1_4, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaPhiDeg4plus, "\\phi _{Fit} - \\phi_{MC}"},
    {kGMTrackDeltaInvQPt, "Global Muon Tracks \\Delta invQPt"},
    {kGMTrackDeltaX, "Global Muon Tracks Delta X at Z_vertex"},
    {kGMTrackDeltaX0_1, "Global Muon Tracks Delta X at Z_vertex"},
    {kGMTrackDeltaX1_4, "Global Muon Tracks Delta X at Z_vertex"},
    {kGMTrackDeltaX4plus, "Global Muon Tracks Delta X at Z_vertex"},
    {kGMTrackDeltaY, "Global Muon Tracks Delta Y at Z_vertex"},
    {kGMTrackR, "Global Muon Tracks Delta R at Z_vertex"},
    {kGMTrackQ, "Global Muon Tracks Charge Match"},
    {kGMTrackQ0_1, "Global Muon Tracks Charge Match (pt < 1)"},
    {kGMTrackQ1_4, "Global Muon Tracks Charge Match (1 < pt < 4)"},
    {kGMTrackQ4plus, "Global Muon Tracks Charge Match (pt > 4)"},
    {kGMTrackChi2, "Global Muon Tracks ~ \\chi^2"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks Pseudorapidity"}};

  std::map<int, std::array<double, 3>> TH1Binning{
    {kGMTracksP, {500, pMin, pMax}},
    {kGMTrackDeltaXErr, {500, -10, 10}},
    {kGMTrackDeltaYErr, {500, -10, 10}},
    {kGMTrackDeltaPhiErr, {500, -10, +10}},
    {kGMTrackDeltaTanLErr, {500, -10, +10}},
    {kGMTrackDeltainvQPtErr, {500, -50, +50}},
    {kMCHResTrackDeltaXErr, {500, -10, 10}},
    {kMCHResTrackDeltaYErr, {500, -10, 10}},
    {kMCHResTrackDeltaPhiErr, {500, -10, +10}},
    {kMCHResTrackDeltaTanLErr, {500, -10, +10}},
    {kMCHResTrackDeltainvQPtErr, {500, -50, +50}},
    {kGMTrackXChi2, {500, 0, 100}},
    {kGMTrackYChi2, {500, 0, 100}},
    {kGMTrackPhiChi2, {500, 0, 100}},
    {kGMTrackTanlChi2, {500, 0, 100}},
    {kGMTrackinvQPtChi2, {500, 0, 100}},
    {kFitChi2, {500, 0, 50}},
    {kGMTrackDeltaTanl, {1000, deltatanlMin, deltatanlMax}},
    {kGMTrackDeltaTanl0_1, {1000, deltatanlMin, deltatanlMax}},
    {kGMTrackDeltaTanl1_4, {1000, deltatanlMin, deltatanlMax}},
    {kGMTrackDeltaTanl4plus, {1000, deltatanlMin, deltatanlMax}},
    {kGMTrackDeltaPhi, {1000, deltaphiMin, deltaphiMax}},
    {kGMTrackDeltaPhi0_1, {1000, deltaphiMin, deltaphiMax}},
    {kGMTrackDeltaPhi1_4, {1000, deltaphiMin, deltaphiMax}},
    {kGMTrackDeltaPhi4plus, {1000, deltaphiMin, deltaphiMax}},
    {kGMTrackDeltaPhiDeg,
     {1000, TMath::RadToDeg() * deltaphiMin,
      TMath::RadToDeg() * deltaphiMax}},
    {kGMTrackDeltaPhiDeg0_1,
     {1000, TMath::RadToDeg() * deltaphiMin,
      TMath::RadToDeg() * deltaphiMax}},
    {kGMTrackDeltaPhiDeg1_4,
     {1000, TMath::RadToDeg() * deltaphiMin,
      TMath::RadToDeg() * deltaphiMax}},
    {kGMTrackDeltaPhiDeg4plus,
     {1000, TMath::RadToDeg() * deltaphiMin,
      TMath::RadToDeg() * deltaphiMax}},
    {kGMTrackDeltaInvQPt, {1000, -10., 10.}},
    {kGMTrackDeltaX, {1000, -.5, .5}},
    {kGMTrackDeltaX0_1, {1000, -.5, .5}},
    {kGMTrackDeltaX1_4, {1000, -.5, .5}},
    {kGMTrackDeltaX4plus, {1000, -.5, .5}},
    {kGMTrackDeltaY, {1000, -.5, .5}},
    {kGMTrackR, {250, 0, 0.5}},
    {kGMTrackQ, {5, -2.1, 2.1}},
    {kGMTrackQ0_1, {5, -2.1, 2.1}},
    {kGMTrackQ1_4, {5, -2.1, 2.1}},
    {kGMTrackQ4plus, {5, -2.1, 2.1}},
    {kGMTrackChi2, {10000, 0, 1000}},
    {kMCTrackspT, {5000, 0, 50}},
    {kMCTracksp, {1000, pMin, pMax}},
    {kMCTrackEta, {1000, etaMin, etaMax}}};

  std::map<int, const char*> TH1XaxisTitles{
    {kGMTracksP, "p [GeV]"},
    {kGMTrackDeltaXErr, "\\Delta x  /\\sigma_{x}"},
    {kGMTrackDeltaYErr, "\\Delta y  /\\sigma_{y}"},
    {kGMTrackDeltaPhiErr, "\\Delta \\phi  /\\sigma_{\\phi}"},
    {kGMTrackDeltaTanLErr, "\\Delta tanl /\\sigma_{tanl}"},
    {kGMTrackDeltainvQPtErr, "\\Delta (q/p_t)/\\sigma_{q/Pt}"},
    {kMCHResTrackDeltaXErr, "\\Delta x  /\\sigma_{x}"},
    {kMCHResTrackDeltaYErr, "\\Delta y  /\\sigma_{y}"},
    {kMCHResTrackDeltaPhiErr, "\\Delta \\phi  /\\sigma_{\\phi}"},
    {kMCHResTrackDeltaTanLErr, "\\Delta tanl /\\sigma_{tanl}"},
    {kMCHResTrackDeltainvQPtErr, "\\Delta (q/p_t)/\\sigma_{q/Pt}"},
    {kGMTrackDeltaTanl, "\\Delta tanl"},
    {kGMTrackXChi2, "\\chi^2"},
    {kGMTrackYChi2, "\\chi^2"},
    {kGMTrackPhiChi2, "\\chi^2"},
    {kGMTrackTanlChi2, "\\chi^2"},
    {kGMTrackinvQPtChi2, "\\chi^2"},
    {kFitChi2, "\\chi^2"},
    {kGMTrackDeltaTanl0_1, "\\Delta tanl"},
    {kGMTrackDeltaTanl1_4, "\\Delta tanl"},
    {kGMTrackDeltaTanl4plus, "\\Delta tanl"},
    {kGMTrackDeltaPhi, "\\Delta \\phi ~[rad]"},
    {kGMTrackDeltaPhi0_1, "\\Delta \\phi ~[rad]"},
    {kGMTrackDeltaPhi1_4, "\\Delta \\phi ~[rad]"},
    {kGMTrackDeltaPhi4plus, "\\Delta \\phi ~[rad]"},
    {kGMTrackDeltaPhiDeg, "\\Delta \\phi ~[deg]"},
    {kGMTrackDeltaPhiDeg0_1, "\\Delta \\phi ~[deg]"},
    {kGMTrackDeltaPhiDeg1_4, "\\Delta \\phi ~[deg]"},
    {kGMTrackDeltaPhiDeg4plus, "\\Delta \\phi ~[deg]"},
    {kGMTrackDeltaInvQPt, "\\Delta invQPt"},
    {kGMTrackDeltaX, "\\Delta x ~[cm]"},
    {kGMTrackDeltaX0_1, "\\Delta x ~[cm]"},
    {kGMTrackDeltaX1_4, "\\Delta x ~[cm]"},
    {kGMTrackDeltaX4plus, "\\Delta x ~[cm]"},
    {kGMTrackDeltaY, "\\Delta y ~[cm]"},
    {kGMTrackR, "\\Delta r ~[cm]"},
    {kGMTrackQ, "q_{fit}-q_{MC}"},
    {kGMTrackQ0_1, "q_{fit}-q_{MC}"},
    {kGMTrackQ1_4, "q_{fit}-q_{MC}"},
    {kGMTrackQ4plus, "q_{fit}-q_{MC}"},
    {kGMTrackChi2, "\\chi^2"},
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
    // h->GetXaxis()->SetLabelSize(0.05);
    // h->GetXaxis()->SetTitleSize(0.05);
    // h->GetYaxis()->SetLabelSize(0.06);
    // h->GetYaxis()->SetTitleSize(0.06);
    h->SetOption("COLZ");
    ++n2Histo;
  }

  // Profiles histograms
  auto PtRes_Profile = new TProfile("Pt_res_prof", "Profile of pt{fit}/pt{MC}",
                                    14, 0, 7, 0, 20, "s");
  PtRes_Profile->GetXaxis()->SetTitle("pt_{MC}");
  PtRes_Profile->GetYaxis()->SetTitle("mean(Pt_{Fit}/Pt_{MC})");

  auto DeltaX_Profile = new TProfile("DeltaX_prof", "Vertexing resolution", 14,
                                     0, 7, -10000., 10000., "s");
  DeltaX_Profile->GetXaxis()->SetTitle("pt_{MC} [GeV]");
  DeltaX_Profile->GetYaxis()->SetTitle("\\sigma_x ~[\\mu m]");

  // TEfficiency histogram
  TEfficiency* qMatchEff = new TEfficiency(
    "QMatchEff", "Charge Match;p_t [GeV];#epsilon", 20, 0, 10);
  // qMatchEff->GetPaintedHistogram()->GetXaxis()->SetLabelSize(0.06);
  // qMatchEff->GetPaintedHistogram()->GetYaxis()->SetLabelSize(0.06);
  // qMatchEff->GetPaintedHistogram()->GetXaxis()->SetTitleSize(0.06);
  // qMatchEff->GetPaintedHistogram()->GetYaxis()->SetTitleSize(0.06);

  TEfficiency* pairedMCHTracksEff = new TEfficiency(
    "PairingEff", "Paired_tracks;p_t [GeV];#epsilon", 20, 0, 10);
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
    "Close_Match_Eff", "Close Matches;p_t [GeV];#epsilon", 20, 0, 10);

  // Counters
  Int_t nChargeMatch = 0;
  Int_t nChargeMiss = 0;
  Int_t nChargeMatch0_1 = 0;
  Int_t nChargeMiss0_1 = 0;
  Int_t nChargeMatch1_4 = 0;
  Int_t nChargeMiss1_4 = 0;
  Int_t nChargeMatch4plus = 0;
  Int_t nChargeMiss4plus = 0;
  Int_t nCorrectMatchGMTracks = 0;
  Int_t nFakeGMTracks = 0;
  Int_t nNoMatchGMTracks = 0;

  // Files & Trees
  // MC
  TFile* o2sim_KineFileIn = new TFile(o2sim_KineFile.c_str());
  TTree* o2SimKineTree = (TTree*)o2sim_KineFileIn->Get("o2sim");

  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  Int_t numberOfEvents = o2SimKineTree->GetEntries();

  // Global Muon Tracks
  TFile* trkFileIn = new TFile(trkFile.c_str());
  TTree* gmTrackTree = (TTree*)trkFileIn->Get("o2sim");
  std::vector<GlobalMuonTrack> trackGMVec, *trackGMVecP = &trackGMVec;
  gmTrackTree->SetBranchAddress("GlobalMuonTrack", &trackGMVecP);

  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mcLabels = nullptr;
  gmTrackTree->SetBranchAddress("GlobalMuonTrackMCTruth", &mcLabels);

  MatchingHelper *matching_helperPtr, matching_helper;
  gDirectory->GetObject("Matching Helper", matching_helperPtr);
  matching_helper = *matching_helperPtr;

  std::string annotation = matching_helper.Annotation();
  std::cout << "matching_helper.Generator = " << matching_helper.Generator
            << std::endl;
  std::cout << "matching_helper.GeneratorConfig = "
            << matching_helper.GeneratorConfig << std::endl;
  std::cout << "matching_helper.MatchingFunction = "
            << matching_helper.MatchingFunction << std::endl;
  std::cout << "matching_helper.MatchingCutFunc = "
            << matching_helper.MatchingCutFunc << std::endl;
  std::cout << "matching_helper.MatchingCutConfig = "
            << matching_helper.MatchingCutConfig << std::endl;
  std::cout << "Annotation = " << annotation << std::endl;

  // MFT Tracks
  TFile* mfttrkFileIn = new TFile("mfttracks.root");
  TTree* mftTrackTree = (TTree*)mfttrkFileIn->Get("o2sim");
  // std::vector<o2::mft::TrackMFT> trackMFTVec, *trackMFTVecP = &trackMFTVec;
  // mftTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);

  std::vector<o2::MCCompLabel>* mftMcLabels = nullptr;
  mftTrackTree->SetBranchAddress("MFTTrackMCTruth", &mftMcLabels);
  mftTrackTree->GetEntry(0);

  gmTrackTree->GetEntry(0);
  o2SimKineTree->GetEntry(0);

  auto field_z = getZField(0, 0, -61.4); // Get field at Center of MFT

  std::string outfilename = "GlobalMuonChecks.root";
  TFile outFile(outfilename.c_str(), "RECREATE");

  // Reconstructed Global Muon Tracks
  std::cout << "Loop over events and reconstructed Global Muon Tracks!"
            << std::endl;
  // GMTracks - Identify reconstructed tracks
  auto nCloseMatches = 0;
  for (int iEvent = 0; iEvent < numberOfEvents; iEvent++) {
    auto iTrack = 0;
    if (DEBUG_VERBOSE) {
      std::cout << "Event = " << iEvent << " with " << trackGMVec.size()
                << " MCH tracks " << std::endl;
    }
    o2SimKineTree->GetEntry(iEvent);
    gmTrackTree->GetEntry(iEvent);

    if (0)
      for (auto& gmTrack : trackGMVec) {
        const auto& label = mcLabels->getLabels(iTrack);
        std::cout << "iTrack = " << iTrack;
        label[0].print();
        iTrack++;
      }

    for (auto& gmTrack : trackGMVec) {

      auto bestMFTTrackMatchID = gmTrack.getBestMFTTrackMatchID();

      const auto& label = mcLabels->getLabels(iTrack);
      // std::cout << "iTrack = " << iTrack;
      // label[0].print();

      if (iEvent == label[0].getEventID()) {
        if (DEBUG_VERBOSE) {
          // std::cout << "  Global Track ID = " <<  iTrack << " ; MFTMatchID =
          // " << bestMFTTrackMatchID << " SourceID = " <<
          // label[0].getSourceID()
          // << " ; EventID = " << label[0].getEventID() << ":  trackID = " <<
          // label[0].getTrackID() << " ; isFake = " << label[0].isFake() << "
          // Label: ";
          std::cout << "  Global Track ID = " << iTrack
                    << " ; MFTMatchID = " << bestMFTTrackMatchID << " Label: ";
          label[0].print();

          // std::cout << "        bestMFTTrackMatchID = " <<
          // bestMFTTrackMatchID << " / labelMFTBestMatch = ";
          // labelMFTBestMatch[0].print();
        }
        if (gmTrack.closeMatch())
          nCloseMatches++;
        pairedMCHTracksEff->Fill(bestMFTTrackMatchID > -1, gmTrack.getPt());
        globalMuonCombinedEff->Fill(label[0].isCorrect(), gmTrack.getPt());
        closeMatchEff->Fill(gmTrack.closeMatch(), gmTrack.getPt());

        if (bestMFTTrackMatchID >= 0) {
          globalMuonCorrectMatchRatio->Fill(label[0].isCorrect(),
                                            gmTrack.getPt());
        }
        if (label[0].isCorrect()) { // Correct match track: add to histograms
          nCorrectMatchGMTracks++;
          // pairedMCHTracksEff->Fill(1,gmTrack.getPt());
          auto thisTrkID = label[0].getTrackID();
          MCTrackT<float>* thisTrack = &(*mcTr).at(thisTrkID);
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
          }

          else {
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
          auto xChi2 = dx * dx / gmTrack.getCovariances()(0, 0);
          auto yChi2 = dy * dy / gmTrack.getCovariances()(1, 1);
          auto phiChi2 = d_Phi * d_Phi / gmTrack.getCovariances()(2, 2);
          auto tanlChi2 = d_tanl * d_tanl / gmTrack.getCovariances()(3, 3);
          auto invQPtChi2 =
            d_invQPt * d_invQPt / sqrt(gmTrack.getCovariances()(4, 4));
          auto fitChi2 = xChi2 + yChi2 + phiChi2 + tanlChi2; // + invQPtChi2;
          auto trackChi2 = gmTrack.getTrackChi2();
          TH1Histos[kGMTracksP]->Fill(gmTrack.getP());
          TH1Histos[kGMTrackDeltaTanl]->Fill(d_tanl);
          TH1Histos[kGMTrackDeltaPhi]->Fill(d_Phi);
          TH1Histos[kGMTrackDeltaInvQPt]->Fill(d_invQPt);
          TH1Histos[kGMTrackDeltaPhiDeg]->Fill(TMath::RadToDeg() * d_Phi);
          TH1Histos[kGMTrackDeltaX]->Fill(dx);

          // std::cout << "DeltaX / sigmaX = " <<
          // dx/sqrt(gmTrack.getCovariances()(0,0)) << std::endl;
          TH1Histos[kGMTrackDeltaXErr]->Fill(
            dx / sqrt(gmTrack.getCovariances()(0, 0)));
          // std::cout << "DeltaY / sigmaY = " <<
          // dy/sqrt(gmTrack.getCovariances()(1,1)) << std::endl;
          TH1Histos[kGMTrackDeltaYErr]->Fill(
            dy / sqrt(gmTrack.getCovariances()(1, 1)));
          // std::cout << "DeltaPhi / sigmaPhi = " <<
          // d_Phi/sqrt(gmTrack.getCovariances()(2,2)) << std::endl;
          TH1Histos[kGMTrackDeltaPhiErr]->Fill(
            d_Phi / sqrt(gmTrack.getCovariances()(2, 2)));
          // std::cout << "DeltaTanl / sigmaTanl = " <<
          // d_tanl/sqrt(gmTrack.getCovariances()(3,3)) << std::endl;
          TH1Histos[kGMTrackDeltaTanLErr]->Fill(
            d_tanl / sqrt(gmTrack.getCovariances()(3, 3)));
          // std::cout << "DeltaPt / sigmaPt = " <<
          // d_Pt/sqrt(gmTrack.getCovariances()(4,4)) << std::endl;
          TH1Histos[kGMTrackDeltainvQPtErr]->Fill(
            d_invQPt / sqrt(gmTrack.getCovariances()(4, 4)));

          //
          TH1Histos[kMCHResTrackDeltaXErr]->Fill(gmTrack.getResiduals2Cov()(0));
          TH1Histos[kMCHResTrackDeltaYErr]->Fill(gmTrack.getResiduals2Cov()(1));
          TH1Histos[kMCHResTrackDeltaPhiErr]->Fill(
            gmTrack.getResiduals2Cov()(2));
          TH1Histos[kMCHResTrackDeltaTanLErr]->Fill(
            gmTrack.getResiduals2Cov()(3));
          TH1Histos[kMCHResTrackDeltainvQPtErr]->Fill(
            gmTrack.getResiduals2Cov()(4));

          TH1Histos[kGMTrackXChi2]->Fill(xChi2);
          TH1Histos[kGMTrackYChi2]->Fill(yChi2);
          TH1Histos[kGMTrackPhiChi2]->Fill(phiChi2);
          TH1Histos[kGMTrackTanlChi2]->Fill(tanlChi2);
          TH1Histos[kGMTrackinvQPtChi2]->Fill(invQPtChi2);
          TH1Histos[kFitChi2]->Fill(fitChi2);
          TH2Histos[kGMTrackChi2vsFitChi2]->Fill(fitChi2, trackChi2);

          DeltaX_Profile->Fill(Pt_MC, dx * 1e4);
          TH1Histos[kGMTrackDeltaY]->Fill(dy);
          TH1Histos[kGMTrackR]->Fill(sqrt(dx * dx + dy * dy));
          TH1Histos[kGMTrackQ]->Fill(d_Charge);
          TH1Histos[kGMTrackChi2]->Fill(trackChi2);
          TH2Histos[kGMTrackDeltaXYVertex]->Fill(10. * dx, 10. * dy);
          TH2Histos[kGMTrackQPRec_MC]->Fill(P_MC * Q_MC, P_fit * Q_fit);
          TH2Histos[kGMTrackPtResolution]->Fill(Pt_MC, Pt_fit / Pt_MC);
          PtRes_Profile->Fill(Pt_MC, Pt_fit / Pt_MC);
          TH2Histos[kGMTrackInvPtResolution]->Fill(
            Pt_MC, (1.0 / Pt_fit - 1.0 / Pt_MC) * Pt_MC);

          // MC histos
          TH1Histos[kMCTrackspT]->Fill(Pt_MC);
          TH1Histos[kMCTracksp]->Fill(P_MC);
          TH1Histos[kMCTrackEta]->Fill(eta_MC);
          TH2Histos[kMCTracksEtaZ]->Fill(vz_MC, eta_MC);

          // Differential histos
          if (Pt_MC <= 1.0) {
            TH2Histos[kGMTrackDeltaXYVertex0_1]->Fill(10. * dx, 10. * dy);
            TH1Histos[kGMTrackDeltaTanl0_1]->Fill(d_tanl);
            TH1Histos[kGMTrackDeltaPhi0_1]->Fill(d_Phi);
            TH1Histos[kGMTrackDeltaPhiDeg0_1]->Fill(TMath::RadToDeg() * d_Phi);
            TH1Histos[kGMTrackDeltaX0_1]->Fill(dx);
            TH1Histos[kGMTrackQ0_1]->Fill(d_Charge);
            d_Charge ? nChargeMiss0_1++ : nChargeMatch0_1++;
          }
          if (Pt_MC > 1.0 and Pt_MC <= 4) {
            TH2Histos[kGMTrackDeltaXYVertex1_4]->Fill(10. * dx, 10. * dy);
            TH1Histos[kGMTrackDeltaTanl1_4]->Fill(d_tanl);
            TH1Histos[kGMTrackDeltaPhi1_4]->Fill(d_Phi);
            TH1Histos[kGMTrackDeltaPhiDeg1_4]->Fill(TMath::RadToDeg() * d_Phi);
            TH1Histos[kGMTrackDeltaX1_4]->Fill(dx);
            TH1Histos[kGMTrackQ1_4]->Fill(d_Charge);
            d_Charge ? nChargeMiss1_4++ : nChargeMatch1_4++;
          }
          if (Pt_MC > 4.0) {
            TH2Histos[kGMTrackDeltaXYVertex4plus]->Fill(10. * dx, 10. * dy);
            TH1Histos[kGMTrackDeltaTanl4plus]->Fill(d_tanl);
            TH1Histos[kGMTrackDeltaPhi4plus]->Fill(d_Phi);
            TH1Histos[kGMTrackDeltaPhiDeg4plus]->Fill(TMath::RadToDeg() *
                                                      d_Phi);
            TH1Histos[kGMTrackDeltaX4plus]->Fill(dx);
            TH1Histos[kGMTrackQ4plus]->Fill(d_Charge);
            d_Charge ? nChargeMiss4plus++ : nChargeMatch4plus++;
          }

          d_Charge ? nChargeMiss++ : nChargeMatch++;
          qMatchEff->Fill(!d_Charge, Pt_MC);
        } else {
          if (bestMFTTrackMatchID >= 0) {
            nFakeGMTracks++;
          } else
            nNoMatchGMTracks++;
        }
      }
      iTrack++;

    } // Loop on GMTracks
  }   // Loop over events

  Int_t nRecoGMTracks = nCorrectMatchGMTracks + nFakeGMTracks;
  Int_t nMCHTracks = nRecoGMTracks + nNoMatchGMTracks;

  // Customize histograms
  TH1Histos[kGMTrackQ]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch,
         100. * nChargeMatch / (nChargeMiss + nChargeMatch)));
  TH1Histos[kGMTrackQ0_1]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch0_1,
         100. * nChargeMatch0_1 / (nChargeMiss0_1 + nChargeMatch0_1)));
  TH1Histos[kGMTrackQ1_4]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch1_4,
         100. * nChargeMatch1_4 / (nChargeMiss1_4 + nChargeMatch1_4)));
  TH1Histos[kGMTrackQ4plus]->SetTitle(
    Form("nChargeMatch = %d (%.2f%%)", nChargeMatch4plus,
         100. * nChargeMatch4plus / (nChargeMiss4plus + nChargeMatch4plus)));

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
  TH2Histos[kGMTrackQPRec_MC]->SetStats(0);
  TH2Histos[kGMTrackPtResolution]->SetStats(0);
  TH2Histos[kGMTrackInvPtResolution]->SetStats(0);
  TH2Histos[kMCTracksEtaZ]->SetStats(0);
  PtRes_Profile->SetStats(0);
  DeltaX_Profile->SetStats(0);
  TH1Histos[kGMTrackQ]->SetStats(0);

  // Fit Slices: Pt resolution
  FitSlicesy(*TH2Histos[kGMTrackInvPtResolution], *TH2Histos[kGMTrackQPRec_MC]);
  FitSlicesy(*TH2Histos[kGMTrackPtResolution], *TH2Histos[kGMTrackQPRec_MC]);

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
    "Matching Summary", annotation, 0, 0, 0, 0, 0, 0, "-", "-", "-",
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex]->GetEntries()),
    "-", "-");

  // Parameters resolution
  auto param_resolution = summary_report_3x2(
    *TH2Histos[kGMTrackDeltaXYVertex], *TH2Histos[kGMTrackPtResolution],
    *PtRes_Profile, *DeltaX_Error, *TH2Histos[kGMTrackQPRec_MC], *qMatchEff,
    "Param Summary", annotation, 0, 0, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackPtResolution]->Integral() /
                     TH2Histos[kGMTrackPtResolution]->GetEntries()),
    "-", "-",
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPRec_MC]->GetEntries()),
    "-");

  // Covariances summary
  auto covariances_summary = summary_report_3x2(
    *TH1Histos[kGMTrackDeltaXErr], *TH1Histos[kGMTrackDeltaPhiErr],
    *TH1Histos[kGMTrackDeltainvQPtErr], *TH1Histos[kGMTrackDeltaYErr],
    *TH1Histos[kGMTrackDeltaTanLErr], *TH2Histos[kGMTrackQPRec_MC],
    "Covariances Summary", annotation, 1, 1, 1, 1, 1, 0,
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaXErr]->Integral() /
                     TH1Histos[kGMTrackDeltaXErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiErr]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltainvQPtErr]->Integral() /
                     TH1Histos[kGMTrackDeltainvQPtErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaYErr]->Integral() /
                     TH1Histos[kGMTrackDeltaYErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaTanLErr]->Integral() /
                     TH1Histos[kGMTrackDeltaTanLErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPRec_MC]->GetEntries()));

  // MCH Residuals Covariances summary
  auto MCHcovariances_summary = summary_report_3x2(
    *TH1Histos[kMCHResTrackDeltaXErr], *TH1Histos[kMCHResTrackDeltaPhiErr],
    *TH1Histos[kMCHResTrackDeltainvQPtErr], *TH1Histos[kMCHResTrackDeltaYErr],
    *TH1Histos[kMCHResTrackDeltaTanLErr], *TH2Histos[kGMTrackQPRec_MC],
    "MCH residuals Covariances Summary", annotation, 1, 1, 1, 1, 1, 0,
    Form("%.2f%%", 100.0 * TH1Histos[kMCHResTrackDeltaXErr]->Integral() /
                     TH1Histos[kMCHResTrackDeltaXErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kMCHResTrackDeltaPhiErr]->Integral() /
                     TH1Histos[kMCHResTrackDeltaPhiErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kMCHResTrackDeltainvQPtErr]->Integral() /
                     TH1Histos[kMCHResTrackDeltainvQPtErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kMCHResTrackDeltaYErr]->Integral() /
                     TH1Histos[kMCHResTrackDeltaYErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kMCHResTrackDeltaTanLErr]->Integral() /
                     TH1Histos[kMCHResTrackDeltaTanLErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPRec_MC]->GetEntries()));

  // Covariances summary 3x3
  auto par_cov_summary3x3 = summary_report_3x3(
    *TH2Histos[kGMTrackDeltaXYVertex], *TH1Histos[kGMTrackDeltaXErr],
    *TH1Histos[kGMTrackDeltaYErr], *DeltaX_Error,
    *TH2Histos[kGMTrackQPRec_MC], *TH1Histos[kGMTrackDeltaPhiErr], *qMatchEff,
    *TH1Histos[kGMTrackDeltainvQPtErr], *TH1Histos[kGMTrackDeltaTanLErr],
    "par_cov_summary3x3", annotation, 0, 1, 1, 0, 0, 1, 0, 1, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaXErr]->Integral() /
                     TH1Histos[kGMTrackDeltaXErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaYErr]->Integral() /
                     TH1Histos[kGMTrackDeltaYErr]->GetEntries()),
    "-",
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPRec_MC]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiErr]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiErr]->GetEntries()),
    "-",
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltainvQPtErr]->Integral() /
                     TH1Histos[kGMTrackDeltainvQPtErr]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaTanLErr]->Integral() /
                     TH1Histos[kGMTrackDeltaTanLErr]->GetEntries()));

  auto param_summary_diff_pt = summary_report_3x3(
    *TH1Histos[kGMTrackDeltaX0_1], *TH1Histos[kGMTrackDeltaTanl0_1],
    *TH1Histos[kGMTrackDeltaPhiDeg0_1], *TH1Histos[kGMTrackDeltaX1_4],
    *TH1Histos[kGMTrackDeltaTanl1_4], *TH1Histos[kGMTrackDeltaPhiDeg1_4],
    *TH1Histos[kGMTrackDeltaX4plus], *TH1Histos[kGMTrackDeltaTanl4plus],
    *TH1Histos[kGMTrackDeltaPhiDeg4plus], "ParamSummaryVsPt", annotation, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaX0_1]->Integral() /
                     TH1Histos[kGMTrackDeltaX0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaTanl0_1]->Integral() /
                     TH1Histos[kGMTrackDeltaTanl0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiDeg0_1]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiDeg0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaX1_4]->Integral() /
                     TH1Histos[kGMTrackDeltaX1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaTanl1_4]->Integral() /
                     TH1Histos[kGMTrackDeltaTanl1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiDeg1_4]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiDeg1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaX4plus]->Integral() /
                     TH1Histos[kGMTrackDeltaX4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaTanl4plus]->Integral() /
                     TH1Histos[kGMTrackDeltaTanl4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiDeg4plus]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiDeg4plus]->GetEntries()));

  auto pt_resolution = summary_report(
    *TH2Histos[kGMTrackPtResolution], *TH2Histos[kGMTrackQPRec_MC],
    *PtRes_Profile, *qMatchEff, "Pt Summary", annotation, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackPtResolution]->Integral() /
                     TH2Histos[kGMTrackPtResolution]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPRec_MC]->GetEntries()));

  auto invpt_resolution = summary_report(
    *TH2Histos[kGMTrackInvPtResolution], *TH2Histos[kGMTrackQPRec_MC],
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kGMTrackInvPtResolution]->GetName()) +
       std::string("_1"))
        .c_str()),
    *(TH1F*)gDirectory->Get(
      (std::string(TH2Histos[kGMTrackInvPtResolution]->GetName()) +
       std::string("_2"))
        .c_str()),
    "InvPt Summary", annotation, 0, 0, 0, 0,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackInvPtResolution]->Integral() /
                     TH2Histos[kGMTrackInvPtResolution]->GetEntries()),
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackQPRec_MC]->Integral() /
                     TH2Histos[kGMTrackQPRec_MC]->GetEntries()));

  auto vertexing_resolution = summary_report(
    *TH2Histos[kGMTrackDeltaXYVertex], *TH1Histos[kGMTrackDeltaX],
    *DeltaX_Error, *TH1Histos[kGMTrackDeltaPhiDeg], "Vertexing Summary",
    annotation, 0, 1, 0, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaX]->Integral() /
                     TH1Histos[kGMTrackDeltaX]->GetEntries()),
    Form("-"),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiDeg]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiDeg]->GetEntries()));

  auto vertexing_resolution0_1 = summary_report(
    *TH2Histos[kGMTrackDeltaXYVertex0_1], *TH1Histos[kGMTrackDeltaX0_1],
    *TH1Histos[kGMTrackDeltaTanl0_1], *TH1Histos[kGMTrackDeltaPhiDeg0_1],
    "Vertexing Summary pt < 1", annotation, 0, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex0_1]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaX0_1]->Integral() /
                     TH1Histos[kGMTrackDeltaX0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaTanl0_1]->Integral() /
                     TH1Histos[kGMTrackDeltaTanl0_1]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiDeg0_1]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiDeg0_1]->GetEntries()));

  auto vertexing_resolution1_4 = summary_report(
    *TH2Histos[kGMTrackDeltaXYVertex1_4], *TH1Histos[kGMTrackDeltaX1_4],
    *TH1Histos[kGMTrackDeltaTanl1_4], *TH1Histos[kGMTrackDeltaPhiDeg1_4],
    "Vertexing Summary 1 < p_t < 4", annotation, 0, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex1_4]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaX1_4]->Integral() /
                     TH1Histos[kGMTrackDeltaX1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaTanl1_4]->Integral() /
                     TH1Histos[kGMTrackDeltaTanl1_4]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiDeg1_4]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiDeg1_4]->GetEntries()));

  auto vertexing_resolution4plus = summary_report(
    *TH2Histos[kGMTrackDeltaXYVertex4plus], *TH1Histos[kGMTrackDeltaX4plus],
    *TH1Histos[kGMTrackDeltaTanl4plus], *TH1Histos[kGMTrackDeltaPhiDeg4plus],
    "Vertexing Summary p_t > 4", annotation, 0, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH2Histos[kGMTrackDeltaXYVertex4plus]->Integral() /
                     TH2Histos[kGMTrackDeltaXYVertex4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaX4plus]->Integral() /
                     TH1Histos[kGMTrackDeltaX4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaTanl4plus]->Integral() /
                     TH1Histos[kGMTrackDeltaTanl4plus]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackDeltaPhiDeg4plus]->Integral() /
                     TH1Histos[kGMTrackDeltaPhiDeg4plus]->GetEntries()));

  auto chi2_summary = summary_report(
    *TH1Histos[kGMTrackChi2], *TH1Histos[kGMTrackXChi2],
    *TH1Histos[kGMTrackTanlChi2], *TH1Histos[kGMTrackPhiChi2], "Chi2 Summary",
    annotation, 1, 1, 1, 1,
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackChi2]->Integral() /
                     TH1Histos[kGMTrackChi2]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackXChi2]->Integral() /
                     TH1Histos[kGMTrackXChi2]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackTanlChi2]->Integral() /
                     TH1Histos[kGMTrackTanlChi2]->GetEntries()),
    Form("%.2f%%", 100.0 * TH1Histos[kGMTrackPhiChi2]->Integral() /
                     TH1Histos[kGMTrackPhiChi2]->GetEntries()));

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
  outFile.WriteObjectAny(&matching_helper, "MatchingHelper", "Matching Helper");

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
  std::cout << " Tanl_StdDev(pt<1) = "
            << TH1Histos[kGMTrackDeltaTanl0_1]->GetStdDev() << std::endl;
  std::cout << " Tanl_StdDev(1<pt<4) = "
            << TH1Histos[kGMTrackDeltaTanl1_4]->GetStdDev() << std::endl;
  std::cout << " Tanl_StdDev(pt>4) = "
            << TH1Histos[kGMTrackDeltaTanl4plus]->GetStdDev() << std::endl;
  std::cout << " Phi_mean = " << TH1Histos[kGMTrackDeltaPhi]->GetMean()
            << std::endl;
  std::cout << " Phi_StdDev = " << TH1Histos[kGMTrackDeltaPhi]->GetStdDev()
            << std::endl;
  std::cout << " Phi_StdDev(pt<1) = "
            << TH1Histos[kGMTrackDeltaPhi0_1]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDev(1<pt<4) = "
            << TH1Histos[kGMTrackDeltaPhi1_4]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDev(pt>4) = "
            << TH1Histos[kGMTrackDeltaPhi4plus]->GetStdDev() << std::endl;
  std::cout << " Phi_meanDeg = " << TH1Histos[kGMTrackDeltaPhiDeg]->GetMean()
            << std::endl;
  std::cout << " Phi_StdDevDeg = "
            << TH1Histos[kGMTrackDeltaPhiDeg]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(pt<1) = "
            << TH1Histos[kGMTrackDeltaPhiDeg0_1]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(1<pt<4) = "
            << TH1Histos[kGMTrackDeltaPhiDeg1_4]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(pt>4) = "
            << TH1Histos[kGMTrackDeltaPhiDeg4plus]->GetStdDev() << std::endl;
  std::cout << " DeltaX_mean = " << TH1Histos[kGMTrackDeltaX]->GetMean()
            << std::endl;
  std::cout << " DeltaX_StdDev = " << TH1Histos[kGMTrackDeltaX]->GetStdDev()
            << std::endl;
  std::cout << " DeltaX_StdDev(pt<1) = "
            << TH1Histos[kGMTrackDeltaX0_1]->GetStdDev() << std::endl;
  std::cout << " DeltaX_StdDev(1<pt<4) = "
            << TH1Histos[kGMTrackDeltaX1_4]->GetStdDev() << std::endl;
  std::cout << " DeltaX_StdDev(pt>4) = "
            << TH1Histos[kGMTrackDeltaX4plus]->GetStdDev() << std::endl;
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
