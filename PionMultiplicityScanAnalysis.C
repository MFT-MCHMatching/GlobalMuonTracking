#if !defined(__CLING__) || defined(__ROOTCLING__)

#endif

#include "TColor.h"
#include "TGraph.h"
#include "include/GlobalMuonTrack.h"
void loadHelperSummaries(std::vector<std::string> fileList,
                         std::vector<float> &NPionsVec,
                         std::vector<float> &EffVec,
                         std::vector<float> &CorrectMatchRatioVec,
                         std::vector<float> &combEffVec,
                         std::vector<float> &closeMatchVec,
                         std::string &matchingConfig);

int get_NPions(std::string);
int get_NEvents(std::string);

void buildEffCanvas(std::vector<std::string> fileList, std::string effName,
                    TFile &, TCanvas &, std::string, std::string);
std::vector<std::string> getSortedFileListByNPions(std::string fileList);

//_________________________________________________________________________________________________
int PionMultiplicityScanAnalysis(
    std::string fileList = "list_of_GMChecks.txt") {

  gStyle->SetFrameLineWidth(3);
  gStyle->SetLineWidth(2);
  gStyle->SetFrameFillColor(10);
  gStyle->SetLabelSize(0.05, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetTitleSize(0.07, "o");
  gStyle->SetStatW(.28);
  gStyle->SetStatH(.26);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPalette(kSolar);
  gStyle->SetMarkerSize(3);

  std::vector<std::string> fileNamesVector =
      getSortedFileListByNPions(fileList);

  std::vector<float> nPionsVec, effVec, combEffVec, correctMatchRatioVec,
      closeMatchVec;
  std::string matchingConfig;
  loadHelperSummaries(fileNamesVector, nPionsVec, effVec, correctMatchRatioVec,
                      combEffVec, closeMatchVec, matchingConfig);

  THStack *effStack = new THStack("effStack", "Matching Efficiencies");
  THStack *correctMatchRatioStack =
      new THStack("correctMatchRatioStack", "GM Tracks CorrectMatchRatio");
  THStack *combinedEffStack = new THStack(
      "combinedEffStack", "Combined Efficiency (eff*correctMatchRatio)");

  TFile *resultsFile = new TFile(
      ("MatchingPerf_PiScan" + matchingConfig + ".root").c_str(), "recreate");
  resultsFile->mkdir("histos");

  TCanvas *cCorrectMatchRatio = new TCanvas("Correct_Match_Ratio_Vs_Pt",
                                            "Correct_Match_Ratio", 1260, 800);
  gStyle->SetOptTitle(0);
  buildEffCanvas(fileNamesVector, "MoreHistos/Correct_Match_Ratio",
                 *resultsFile, *cCorrectMatchRatio, matchingConfig,
                 "Correct_Match_Ratio");

  TPaveText *pt = new TPaveText(0.1, 0.918, 0.9, 0.995, "NDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(4000);
  pt->AddText("Correct_Match_Ratio");
  pt->AddText(matchingConfig.c_str());
  pt->Draw();
  // gSystem->ProcessEvents();
  // cCorrectMatchRatio->Update();

  TCanvas *cEff =
      new TCanvas("PairingEff_Vs_Pt", "Pairing_efficiency", 1260, 800);
  gStyle->SetOptTitle(0);
  buildEffCanvas(fileNamesVector, "MoreHistos/PairingEff", *resultsFile, *cEff,
                 matchingConfig, "Pairing_Efficiency");
  pt = new TPaveText(0.1, 0.918, 0.9, 0.995, "NDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(4000);
  pt->AddText("Pairing_efficiency");
  pt->AddText(matchingConfig.c_str());
  pt->Draw();
  cEff->Draw();
  cEff->SaveAs(("Pairing_efficiency" + matchingConfig + ".png").c_str());

  TCanvas *cCloseMatch =
      new TCanvas("CloseMatch_Vs_Pt", "Close_Match", 1260, 800);
  gStyle->SetOptTitle(0);
  buildEffCanvas(fileNamesVector, "MoreHistos/Close_Match_Eff", *resultsFile,
                 *cCloseMatch, matchingConfig, "Close_Match");
  pt = new TPaveText(0.1, 0.918, 0.9, 0.995, "NDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(4000);
  pt->AddText("Close_Match");
  pt->AddText(matchingConfig.c_str());
  pt->Draw();
  cCloseMatch->Draw();
  cCloseMatch->SaveAs(("Close_Match" + matchingConfig + ".png").c_str());

  TCanvas *canvasMult_Scan = new TCanvas("Efficiency_Vs_Multiplicity",
                                         "Pion multiplicity scan", 1260, 800);
  TMultiGraph *MultiGraph_Mult_Scan = new TMultiGraph();
  MultiGraph_Mult_Scan->SetName("PiMultScanSummaryMultiGraphs");
  MultiGraph_Mult_Scan->SetTitle(matchingConfig.c_str());

  TGraph *gEff = new TGraph(nPionsVec.size(), &nPionsVec[0], &effVec[0]);
  gEff->SetTitle("Pairing_Efficiency");
  gEff->SetName("PairingEfficiency");
  gEff->GetXaxis()->SetTitle("Pion multiplicity");
  gEff->SetMarkerStyle(20);
  gEff->SetMarkerColor(1);
  gEff->SetMarkerSize(3);

  TGraph *gCorrectMatchRatio =
      new TGraph(nPionsVec.size(), &nPionsVec[0], &correctMatchRatioVec[0]);
  gCorrectMatchRatio->SetTitle("Correct_Match_Ratio");
  gCorrectMatchRatio->SetName("GMTracksCorrectMatchRatio");
  gCorrectMatchRatio->GetXaxis()->SetTitle("Pion multiplicity");
  gCorrectMatchRatio->SetMarkerStyle(21);
  gCorrectMatchRatio->SetMarkerColor(2);
  gCorrectMatchRatio->SetMarkerSize(3);

  TGraph *gCombEff =
      new TGraph(nPionsVec.size(), &nPionsVec[0], &combEffVec[0]);
  gCombEff->SetTitle("GMTracking_Combined_Efficiency");
  gCombEff->SetName("GMTrackingCombinedEfficiency");
  gCombEff->GetXaxis()->SetTitle("Pion multiplicity");
  gCombEff->SetMarkerStyle(22);
  gCombEff->SetMarkerColor(3);
  gCombEff->SetMarkerSize(3);

  TGraph *gCloseMatch =
      new TGraph(nPionsVec.size(), &nPionsVec[0], &closeMatchVec[0]);
  gCloseMatch->SetTitle("GMTracking_Close_Match");
  gCloseMatch->SetName("GMTrackingCloseMatch");
  gCloseMatch->GetXaxis()->SetTitle("Pion multiplicity");
  gCloseMatch->SetMarkerStyle(23);
  gCloseMatch->SetMarkerColor(4);
  gCloseMatch->SetMarkerSize(3);

  MultiGraph_Mult_Scan->Add(gEff);
  MultiGraph_Mult_Scan->Add(gCorrectMatchRatio);
  MultiGraph_Mult_Scan->Add(gCombEff);
  MultiGraph_Mult_Scan->Add(gCloseMatch);

  MultiGraph_Mult_Scan->GetXaxis()->SetTitle("Pion multiplicity");
  MultiGraph_Mult_Scan->SetTitle("Scan");

  gEff->Draw();
  gPad->Modified();
  MultiGraph_Mult_Scan->Draw("LP same");
  MultiGraph_Mult_Scan->GetYaxis()->SetLimits(0., 1.0);
  MultiGraph_Mult_Scan->SetMinimum(0.);
  MultiGraph_Mult_Scan->SetMaximum(1.);

  canvasMult_Scan->Update();
  pt = new TPaveText(0.1, 0.918, 0.9, 0.995, "NDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(4000);
  pt->AddText("Pion multiplicity scan summary");
  pt->AddText(matchingConfig.c_str());
  pt->Draw();

  resultsFile->cd();

  canvasMult_Scan->Draw();
  canvasMult_Scan->BuildLegend();

  canvasMult_Scan->SaveAs(
      ("PiMultScanSummary" + matchingConfig + ".png").c_str());
  canvasMult_Scan->Write();
  resultsFile->cd("histos");
  MultiGraph_Mult_Scan->Write();
  gEff->Write();
  gCorrectMatchRatio->Write();
  gCombEff->Write();
  cCloseMatch->Write();
  resultsFile->Close();
  return 0;
}

//_____________________________________________________________________________________
void loadHelperSummaries(std::vector<std::string> fileList,
                         std::vector<float> &NPionsVec,
                         std::vector<float> &EffVec,
                         std::vector<float> &CorrectMatchRatioVec,
                         std::vector<float> &combEffVec,
                         std::vector<float> &closeMatchVec,
                         std::string &matchingConfig) {

  TFile *checksFile;

  for (auto checksFileName : fileList) {
    checksFile = new TFile(checksFileName.c_str(), "read");
    if (checksFile) {
      MatchingHelper *matching_helperPtr, matching_helper;
      gDirectory->GetObject("Matching Helper", matching_helperPtr);
      matching_helper = *matching_helperPtr;
      int nPions = get_NPions(checksFileName);
      int nEvents = get_NEvents(checksFileName);
      NPionsVec.push_back(nPions);
      auto correctMatchRatio = matching_helper.getCorrectMatchRatio();
      auto eff = matching_helper.getPairingEfficiency();
      CorrectMatchRatioVec.push_back(correctMatchRatio);
      EffVec.push_back(eff);
      combEffVec.push_back(eff * correctMatchRatio);
      closeMatchVec.push_back((1.0 * matching_helper.nCloseMatches) /
                              matching_helper.nMCHTracks);
      matchingConfig = matching_helper.MatchingConfig();
    }
  }
}

//_________________________________________________________________________________________________
void buildEffCanvas(std::vector<std::string> fileList, std::string effName,
                    TFile &resultsFile, TCanvas &canvas,
                    std::string matchingConfig, std::string title) {

  TFile *checksFile;
  int n = 0;
  for (auto checksFileName : fileList) {
    checksFile = new TFile(checksFileName.c_str(), "read");
    canvas.cd();

    if (checksFile) {
      int nPions = get_NPions(checksFileName);
      int nEvents = get_NEvents(checksFileName);
      checksFile->cd("MoreHistos");
      TEfficiency *pEff = (TEfficiency *)checksFile->Get(effName.c_str());
      if (pEff) {
        gStyle->SetOptTitle(0);
        pEff->SetTitle((std::to_string(nPions) + " pions").c_str());
        pEff->SetName(
            (std::string(pEff->GetName()) + std::to_string(nPions) + "_pions")
                .c_str());
        resultsFile.cd("histos");
        pEff->SetMarkerColor(n + 1);
        pEff->SetMarkerSize(2);
        if (!n) {

          canvas.cd();
          pEff->Draw();
          gPad->Update();
          auto graph = pEff->GetPaintedGraph();
          graph->SetMinimum(0);
          graph->SetMaximum(1);
          gPad->Update();

        } else {
          pEff->Draw("SAME");
        }
        pEff->Write();
        gPad->Update();
        auto graph = pEff->GetPaintedGraph();
        graph->SetMinimum(0);
        graph->SetMaximum(1);
        gPad->Update();
        n++;
        // frame->GetYaxis()->SetLimits(0., 1.);
        gSystem->ProcessEvents();
        canvas.Update();
      }
      checksFile->Close();
      delete checksFile;
    }
  }

  resultsFile.cd();
  canvas.cd();

  TPaveText *pt = new TPaveText(0.1, 0.918, 0.9, 0.995, "blNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(4000);
  pt->SetLineWidth(2);
  pt->SetTextFont(42);
  pt->AddText(title.c_str());
  pt->AddText(matchingConfig.c_str());
  pt->Draw();
  canvas.Draw();
  canvas.BuildLegend();
  canvas.SaveAs((title + matchingConfig + ".png").c_str());
  resultsFile.cd();
  canvas.Write();
}

//_________________________________________________________________________________________________
std::vector<std::string> getSortedFileListByNPions(std::string fileList) {
  std::vector<std::string> fileNamesVector;

  std::ifstream listOfGMChecksFile(fileList.c_str());
  TFile *checksFile;
  int n = 0;
  if (listOfGMChecksFile) {
    string checksFileName;
    while (getline(listOfGMChecksFile, checksFileName)) {
      fileNamesVector.push_back(checksFileName);
    }
  }

  struct {
    bool operator()(std::string a, std::string b) const {
      return get_NPions(a) < get_NPions(b);
    }
  } compareNPions;
  std::sort(fileNamesVector.begin(), fileNamesVector.end(), compareNPions);
  std::cout << "Sorted List of filenames by nPions: " << std::endl;
  for (auto filename : fileNamesVector) {
    std::cout << filename << std::endl;
  }

  std::cout << std::endl;
  return fileNamesVector;
}

//_________________________________________________________________________________________________
int get_NEvents(std::string annotation) {
  std::string N = annotation.substr(0, annotation.find("Ev"));
  N = N.substr(N.find_last_of("_") + 1);
  return std::stoi(N);
}

//_________________________________________________________________________________________________
int get_NPions(std::string annotation) {
  std::string N = annotation.substr(0, annotation.find("Pi_"));
  N = N.substr(N.find_last_of("_") + 1);
  return std::stoi(N);
}
