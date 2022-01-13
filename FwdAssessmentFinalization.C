// Development macro to finalize GlobalFwd Assessment

static constexpr std::array<short, 7> sMinNClustersList = {4, 5, 6, 7, 8, 9, 10};
bool mUseMC = true;
enum mMFTTrackTypes
{
  kReco,
  kGen,
  kPairable,
  kRecoTrue,
  kNumberOfTrackTypes
};

// Histos for reconstructed tracks
std::unique_ptr<TH1F> mTrackNumberOfClusters = nullptr;
std::unique_ptr<TH1F> mTrackInvQPt = nullptr;
std::unique_ptr<TH1F> mTrackChi2 = nullptr;
std::unique_ptr<TH1F> mTrackCharge = nullptr;
std::unique_ptr<TH1F> mTrackPhi = nullptr;
std::unique_ptr<TH1F> mTrackEta = nullptr;
std::array<std::unique_ptr<TH1F>, 7> mTrackEtaNCls = {nullptr};
std::array<std::unique_ptr<TH1F>, 7> mTrackPhiNCls = {nullptr};
std::array<std::unique_ptr<TH2F>, 7> mTrackXYNCls = {nullptr};
std::array<std::unique_ptr<TH2F>, 7> mTrackEtaPhiNCls = {nullptr};
std::unique_ptr<TH1F> mTrackTanl = nullptr;

// Histos and data for MC analysis
std::vector<std::string> mNameOfTrackTypes = {"Rec",
                                              "Gen",
                                              "Pairable",
                                              "RecoTrue"};

std::unique_ptr<TH2F> mHistPhiRecVsPhiGen = nullptr;
std::unique_ptr<TH2F> mHistEtaRecVsEtaGen = nullptr;

std::array<std::unique_ptr<TH2F>, kNumberOfTrackTypes> mHistPhiVsEta;
std::array<std::unique_ptr<TH2F>, kNumberOfTrackTypes> mHistPtVsEta;
std::array<std::unique_ptr<TH2F>, kNumberOfTrackTypes> mHistPhiVsPt;
std::array<std::unique_ptr<TH2F>, kNumberOfTrackTypes> mHistZvtxVsEta;
std::array<std::unique_ptr<TH2F>, kNumberOfTrackTypes> mHistRVsZ;

// Histos for reconstruction assessment

std::unique_ptr<TEfficiency> mChargeMatchEff = nullptr;
std::unique_ptr<TEfficiency> mPurityPt = nullptr;
std::unique_ptr<TEfficiency> mPurityPtInner = nullptr;
std::unique_ptr<TEfficiency> mPurityPtOuter = nullptr;
std::unique_ptr<TH2D> mPairingEtaPt = nullptr;

std::vector<std::unique_ptr<TEfficiency>> mPurityPtInnerVec;
std::vector<std::unique_ptr<TEfficiency>> mPurityPtOuterVec;
std::vector<std::unique_ptr<TH2D>> mPurityPtInnerVecTH2;
std::vector<std::unique_ptr<TH2D>> mPurityPtOuterVecTH2;
std::vector<std::unique_ptr<TH1D>> mPairingPtInnerVecTH1;
std::vector<std::unique_ptr<TH1D>> mPairingPtOuterVecTH1;
std::vector<std::unique_ptr<TH2D>> mPairingEtaPtVec;

enum TH3HistosCodes
{
  kTH3GMTrackDeltaXDeltaYEta,
  kTH3GMTrackDeltaXDeltaYPt,
  kTH3GMTrackDeltaXVertexPtEta,
  kTH3GMTrackDeltaYVertexPtEta,
  kTH3GMTrackInvQPtResolutionPtEta,
  kTH3GMTrackInvQPtResMCHPtEta,
  kTH3GMTrackXPullPtEta,
  kTH3GMTrackYPullPtEta,
  kTH3GMTrackPhiPullPtEta,
  kTH3GMTrackTanlPullPtEta,
  kTH3GMTrackInvQPtPullPtEta,
  kTH3GMTrackReducedChi2PtEta,
  kTH3GMTrackPtEtaChi2,
  kTH3GMTrackPtEtaMatchScore,
  kTH3GMTruePtEtaChi2,
  kTH3GMTruePtEtaMatchScore,
  kTH3GMCloseMatchPtEtaChi2,
  kTH3GMCloseMatchPtEtaMatchScore,
  kTH3GMPairablePtEtaZ,
  kNTH3Histos
};

std::map<int, const char *> TH3Names{
    {kTH3GMTrackDeltaXDeltaYEta, "TH3GMTrackDeltaXDeltaYEta"},
    {kTH3GMTrackDeltaXDeltaYPt, "TH3GMTrackDeltaXDeltaYPt"},
    {kTH3GMTrackDeltaXVertexPtEta, "TH3GMTrackDeltaXVertexPtEta"},
    {kTH3GMTrackDeltaYVertexPtEta, "TH3GMTrackDeltaYVertexPtEta"},
    {kTH3GMTrackInvQPtResolutionPtEta, "TH3GMTrackInvQPtResolutionPtEta"},
    {kTH3GMTrackInvQPtResMCHPtEta, "TH3GMTrackInvQPtResMCHPtEta"},
    {kTH3GMTrackXPullPtEta, "TH3GMTrackXPullPtEta"},
    {kTH3GMTrackYPullPtEta, "TH3GMTrackYPullPtEta"},
    {kTH3GMTrackPhiPullPtEta, "TH3GMTrackPhiPullPtEta"},
    {kTH3GMTrackTanlPullPtEta, "TH3GMTrackTanlPullPtEta"},
    {kTH3GMTrackInvQPtPullPtEta, "TH3GMTrackInvQPtPullPtEta"},
    {kTH3GMTrackReducedChi2PtEta, "TH3GMTrackReducedChi2PtEta"},
    {kTH3GMCloseMatchPtEtaChi2, "TH3GMCloseMatchPtEtaChi2"},
    {kTH3GMCloseMatchPtEtaMatchScore, "TH3GMCloseMatchPtEtaMatchScore"},
    {kTH3GMPairablePtEtaZ, "TH3GMPairablePtEtaZ"},
    {kTH3GMTrackPtEtaChi2, "TH3GMTrackPtEtaChi2"},
    {kTH3GMTrackPtEtaMatchScore, "TH3GMTrackPtEtaMatchScore"},
    {kTH3GMTruePtEtaChi2, "TH3GMTruePtEtaChi2"},
    {kTH3GMTruePtEtaMatchScore, "TH3GMTruePtEtaMatchScore"}};

std::map<int, const char *> TH3Titles{
    {kTH3GMTrackDeltaXDeltaYEta, "TH3GMTrackDeltaXDeltaYEta"},
    {kTH3GMTrackDeltaXDeltaYPt, "TH3GMTrackDeltaXDeltaYPt"},
    {kTH3GMTrackDeltaXVertexPtEta, "TH3GMTrackDeltaXVertexPtEta"},
    {kTH3GMTrackDeltaYVertexPtEta, "TH3GMTrackDeltaYVertexPtEta"},
    {kTH3GMTrackInvQPtResolutionPtEta, "TH3GMTrackInvQPtResolutionPtEta"},
    {kTH3GMTrackInvQPtResMCHPtEta, "TH3GMTrackInvQPtResMCHPtEta"},
    {kTH3GMTrackXPullPtEta, "TH3GMTrackXPullPtEta"},
    {kTH3GMTrackYPullPtEta, "TH3GMTrackYPullPtEta"},
    {kTH3GMTrackPhiPullPtEta, "TH3GMTrackPhiPullPtEta"},
    {kTH3GMTrackTanlPullPtEta, "TH3GMTrackTanlPullPtEta"},
    {kTH3GMTrackInvQPtPullPtEta, "TH3GMTrackInvQPtPullPtEta"},
    {kTH3GMTrackReducedChi2PtEta, "TH3GMTrackReducedChi2PtEta"},
    {kTH3GMCloseMatchPtEtaChi2, "TH3GMCloseMatchPtEtaChi2"},
    {kTH3GMCloseMatchPtEtaMatchScore, "TH3GMCloseMatchPtEtaMatchScore"},
    {kTH3GMPairablePtEtaZ, "TH3GMPairablePtEtaZ"},
    {kTH3GMTrackPtEtaChi2, "TH3GMTrackPtEtaChi2"},
    {kTH3GMTrackPtEtaMatchScore, "TH3GMTrackPtEtaMatchScore"},
    {kTH3GMTruePtEtaChi2, "TH3GMTruePtEtaChi2"},
    {kTH3GMTruePtEtaMatchScore, "TH3GMTruePtEtaMatchScore"}};

std::map<int, std::array<double, 9>> TH3Binning{
    {kTH3GMTrackDeltaXDeltaYEta, {16, 2.2, 3.8, 1000, -1000, 1000, 1000, -1000, 1000}},
    {kTH3GMTrackDeltaXDeltaYPt, {40, 0, 20, 1000, -1000, 1000, 1000, -1000, 1000}},
    {kTH3GMTrackDeltaYVertexPtEta, {40, 0, 20, 16, 2.2, 3.8, 1000, -1000, 1000}},
    {kTH3GMTrackDeltaXVertexPtEta, {40, 0, 20, 16, 2.2, 3.8, 1000, -1000, 1000}},
    {kTH3GMTrackInvQPtResolutionPtEta, {40, 0, 20, 16, 2.2, 3.8, 1000, -5, 5}},
    {kTH3GMTrackInvQPtResMCHPtEta, {40, 0, 20, 16, 2.2, 3.8, 1000, -5, 5}},
    {kTH3GMTrackXPullPtEta, {40, 0, 20, 16, 2.2, 3.8, 200, -10, 10}},
    {kTH3GMTrackYPullPtEta, {40, 0, 20, 16, 2.2, 3.8, 200, -10, 10}},
    {kTH3GMTrackPhiPullPtEta, {40, 0, 20, 16, 2.2, 3.8, 200, -10, 10}},
    {kTH3GMTrackTanlPullPtEta, {40, 0, 20, 16, 2.2, 3.8, 200, -10, 10}},
    {kTH3GMTrackInvQPtPullPtEta, {40, 0, 20, 16, 2.2, 3.8, 200, -50, 50}},
    {kTH3GMTrackReducedChi2PtEta, {40, 0, 20, 16, 2.2, 3.8, 1000, 0, 100}},
    {kTH3GMCloseMatchPtEtaChi2, {40, 0, 20, 16, 2.2, 3.8, 1000, 0, 100}},
    {kTH3GMCloseMatchPtEtaMatchScore, {40, 0, 20, 16, 2.2, 3.8, 2000, 0, 20.0}},
    {kTH3GMPairablePtEtaZ, {40, 0, 20, 16, 2.2, 3.8, 30, -15, 15}},
    {kTH3GMTrackPtEtaChi2, {40, 0, 20, 16, 2.2, 3.8, 1000, 0, 100}},
    {kTH3GMTrackPtEtaMatchScore, {40, 0, 20, 16, 2.2, 3.8, 2000, 0, 20.0}},
    {kTH3GMTruePtEtaChi2, {40, 0, 20, 16, 2.2, 3.8, 1000, 0, 100}},
    {kTH3GMTruePtEtaMatchScore, {40, 0, 20, 16, 2.2, 3.8, 2000, 0, 20.0}}};

std::map<int, const char *> TH3XaxisTitles{
    {kTH3GMTrackDeltaXDeltaYEta, R"(\\eta_{MC})"},
    {kTH3GMTrackDeltaXDeltaYPt, R"(p_{t}_{MC})"},
    {kTH3GMTrackDeltaXVertexPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackDeltaYVertexPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackInvQPtResolutionPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackInvQPtResMCHPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackXPullPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackYPullPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackPhiPullPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackTanlPullPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackInvQPtPullPtEta, R"(p_{t}_{MC})"},
    {kTH3GMTrackReducedChi2PtEta, R"(p_{t}_{MC})"},
    {kTH3GMCloseMatchPtEtaChi2, R"(p_{t}_{Fit}_{MC})"},
    {kTH3GMCloseMatchPtEtaMatchScore, R"(p_{t}_{Fit}_{MC})"},
    {kTH3GMPairablePtEtaZ, R"(p_{t}_{Fit}_{MC})"},
    {kTH3GMTrackPtEtaChi2, R"(p_{t}_{Fit}_{MC})"},
    {kTH3GMTrackPtEtaMatchScore, R"(p_{t}_{Fit}_{MC})"},
    {kTH3GMTruePtEtaChi2, R"(p_{t}_{Fit}_{MC})"},
    {kTH3GMTruePtEtaMatchScore, R"(p_{t}_{Fit}_{MC})"}};

std::map<int, const char *> TH3YaxisTitles{
    {kTH3GMTrackDeltaXDeltaYEta, R"(X_{residual \rightarrow vtx} (\mu m))"},
    {kTH3GMTrackDeltaXDeltaYPt, R"(X_{residual \rightarrow vtx} (\mu m))"},
    {kTH3GMTrackDeltaXVertexPtEta, R"(\eta_{MC}v)"},
    {kTH3GMTrackDeltaYVertexPtEta, R"(\eta_{MC})"},
    {kTH3GMTrackInvQPtResolutionPtEta, R"(\eta_{MC})"},
    {kTH3GMTrackInvQPtResMCHPtEta, R"(\eta_{MC})"},
    {kTH3GMTrackXPullPtEta, R"(\eta_{MC})"},
    {kTH3GMTrackYPullPtEta, R"(\eta_{MC})"},
    {kTH3GMTrackPhiPullPtEta, R"(\eta_{MC})"},
    {kTH3GMTrackTanlPullPtEta, R"(\eta_{MC})"},
    {kTH3GMTrackInvQPtPullPtEta, R"(\eta_{MC})"},
    {kTH3GMTrackReducedChi2PtEta, R"(\eta_{MC})"},
    {kTH3GMCloseMatchPtEtaChi2, R"(\eta_{Fit})"},
    {kTH3GMCloseMatchPtEtaMatchScore, R"(\eta_{Fit})"},
    {kTH3GMPairablePtEtaZ, R"(\eta_{MC})"},
    {kTH3GMTrackPtEtaChi2, R"(\eta_{Fit})"},
    {kTH3GMTrackPtEtaMatchScore, R"(\eta_{Fit})"},
    {kTH3GMTruePtEtaChi2, R"(\eta_{Fit})"},
    {kTH3GMTruePtEtaMatchScore, R"(\eta_{Fit})"}};

std::map<int, const char *> TH3ZaxisTitles{
    {kTH3GMTrackDeltaXDeltaYEta, R"(Y_{residual \rightarrow vtx} (\mu m))"},
    {kTH3GMTrackDeltaXDeltaYPt, R"(Y_{residual \rightarrow vtx} (\mu m))"},
    {kTH3GMTrackDeltaXVertexPtEta, R"(X_{residual \rightarrow vtx} (\mu m))"},
    {kTH3GMTrackDeltaYVertexPtEta, R"(Y_{residual \rightarrow vtx} (\mu m))"},
    {kTH3GMTrackInvQPtResolutionPtEta, R"((q/p_{t})_{residual}/(q/p_{t}))"},
    {kTH3GMTrackInvQPtResMCHPtEta, R"((q/p_{t})_{residual}/(q/p_{t}))"},
    {kTH3GMTrackXPullPtEta, R"(\Delta X/\sigma_{X})"},
    {kTH3GMTrackYPullPtEta, R"(\Delta Y/\sigma_{Y})"},
    {kTH3GMTrackPhiPullPtEta, R"(\Delta \phi/\sigma_{\phi})"},
    {kTH3GMTrackTanlPullPtEta, R"(\Delta \tan(\lambda)/\sigma_{tan(\lambda)})"},
    {kTH3GMTrackInvQPtPullPtEta, R"((\Delta q/p_t)/\sigma_{q/p_{t}})"},
    {kTH3GMTrackReducedChi2PtEta, R"(\chi^2/d.f.)"},
    {kTH3GMCloseMatchPtEtaChi2, R"(Match \chi^2)"},
    {kTH3GMCloseMatchPtEtaMatchScore, R"(Matching Score)"},
    {kTH3GMPairablePtEtaZ, R"(z_{vtx})"},
    {kTH3GMTrackPtEtaChi2, R"(Match \chi^2)"},
    {kTH3GMTrackPtEtaMatchScore, R"(Matching Score)"},
    {kTH3GMTruePtEtaChi2, R"(Match \chi^2)"},
    {kTH3GMTruePtEtaMatchScore, R"(Matching Score)"}};

enum TH3SlicedCodes
{
  kDeltaXVertexVsEta,
  kDeltaXVertexVsPt,
  kDeltaYVertexVsEta,
  kDeltaYVertexVsPt,
  kXPullVsEta,
  kXPullVsPt,
  kYPullVsEta,
  kYPullVsPt,
  kInvQPtResVsEta,
  kInvQPtResVsPt,
  kInvQPtResMCHVsEta,
  kInvQPtResMCHVsPt,
  kPhiPullVsEta,
  kPhiPullVsPt,
  kTanlPullVsEta,
  kTanlPullVsPt,
  kInvQPtPullVsEta,
  kInvQPtPullVsPt,
  kNSlicedTH3
};

std::map<int, const char *> TH3SlicedNames{
    {kDeltaXVertexVsEta, "DeltaXVertexVsEta"},
    {kDeltaXVertexVsPt, "DeltaXVertexVsPt"},
    {kDeltaYVertexVsEta, "DeltaYVertexVsEta"},
    {kDeltaYVertexVsPt, "DeltaYVertexVsPt"},
    {kXPullVsEta, "XPullVsEta"},
    {kXPullVsPt, "XPullVsPt"},
    {kYPullVsEta, "YPullVsEta"},
    {kYPullVsPt, "YPullVsPt"},
    {kInvQPtResVsEta, "InvQPtResVsEta"},
    {kInvQPtResVsPt, "InvQPtResVsPt"},
    {kInvQPtResMCHVsEta, "InvQPtResMCHVsEta"},
    {kInvQPtResMCHVsPt, "InvQPtResMCHVsPt"},
    {kPhiPullVsEta, "PhiPullVsEta"},
    {kPhiPullVsPt, "PhiPullVsPt"},
    {kTanlPullVsEta, "TanlPullVsEta"},
    {kTanlPullVsPt, "TanlPullVsPt"},
    {kInvQPtPullVsEta, "InvQPtPullVsEta"},
    {kInvQPtPullVsPt, "InvQPtPullVsPt"}};

std::map<int, int> TH3SlicedMap{
    {kDeltaXVertexVsEta, kTH3GMTrackDeltaXVertexPtEta},
    {kDeltaXVertexVsPt, kTH3GMTrackDeltaXVertexPtEta},
    {kDeltaYVertexVsEta, kTH3GMTrackDeltaYVertexPtEta},
    {kDeltaYVertexVsPt, kTH3GMTrackDeltaYVertexPtEta},
    {kXPullVsEta, kTH3GMTrackXPullPtEta},
    {kXPullVsPt, kTH3GMTrackXPullPtEta},
    {kYPullVsEta, kTH3GMTrackYPullPtEta},
    {kYPullVsPt, kTH3GMTrackYPullPtEta},
    {kInvQPtResVsEta, kTH3GMTrackInvQPtResolutionPtEta},
    {kInvQPtResVsPt, kTH3GMTrackInvQPtResolutionPtEta},
    {kInvQPtResMCHVsEta, kTH3GMTrackInvQPtResMCHPtEta},
    {kInvQPtResMCHVsPt, kTH3GMTrackInvQPtResMCHPtEta},
    {kPhiPullVsEta, kTH3GMTrackPhiPullPtEta},
    {kPhiPullVsPt, kTH3GMTrackPhiPullPtEta},
    {kTanlPullVsEta, kTH3GMTrackTanlPullPtEta},
    {kTanlPullVsPt, kTH3GMTrackTanlPullPtEta},
    {kInvQPtPullVsEta, kTH3GMTrackInvQPtPullPtEta},
    {kInvQPtPullVsPt, kTH3GMTrackInvQPtPullPtEta}};

std::array<std::unique_ptr<TH3F>, kNTH3Histos> mTH3Histos;
std::array<TCanvas *, kNSlicedTH3> mSlicedCanvas;

enum GMAssesmentCanvases
{
  kPurityPtOuter,
  kPurityPtInner,
  kPairingEffPtOuter,
  kPairingEffPtInner,
  kPurityVsEfficiencyVeryLowPt,
  kPurityVsEfficiencyIntegrated,
  kNGMAssesmentCanvases
};

std::map<int, const char *> GMAssesmentNames{
    {kPurityPtOuter, "PurityPtOuter"},
    {kPurityPtInner, "PurityPtInner"},
    {kPairingEffPtOuter, "PairingEffPtOuter"},
    {kPairingEffPtInner, "PairingEffPtInner"},
    {kPurityVsEfficiencyVeryLowPt, "PurityVsEfficiencyVeryLowPt"},
    {kPurityVsEfficiencyIntegrated, "PurityVsEfficiencyIntegrated"}};

std::array<TCanvas *, kNGMAssesmentCanvases> mAssessmentCanvas;
void TH3Slicer(TCanvas *canvas, std::unique_ptr<TH3F> &histo3D, std::vector<float> list, double window, int iPar, float marker_size = 1.5);

//-------------------------------------------------------------------------------------------
void TH3Slicer(TCanvas *canvas, std::unique_ptr<TH3F> &histo3D, std::vector<float> list, double window, int iPar, float marker_size)
{
  std::string cname = canvas->GetName();
  std::string ctitle = cname;
  std::string option;

  TObjArray aSlices;
  histo3D->GetYaxis()->SetRange(0, 0);
  histo3D->GetXaxis()->SetRange(0, 0);
  bool first = true;
  if (cname.find("VsEta") < cname.length())
  {
    for (auto ptmin : list)
    {
      auto ptmax = ptmin + window;
      histo3D->GetXaxis()->SetRangeUser(ptmin, ptmax);

      std::string ytitle = "\\sigma (";
      ytitle += histo3D->GetZaxis()->GetTitle();
      ytitle += ")";
      auto title = Form("_%1.2f_%1.2f_yz", ptmin, ptmax);
      auto aDBG = (TH2F *)histo3D->Project3D(title);
      aDBG->GetXaxis()->SetRangeUser(0, 0);

      aDBG->FitSlicesX(nullptr, 0, -1, 4, "QNR", &aSlices);
      auto th1DBG = (TH1F *)aSlices[iPar];
      th1DBG->SetTitle(Form("%1.2f < p_t < %1.2f", ptmin, ptmax));
      th1DBG->SetStats(0);
      th1DBG->SetYTitle(ytitle.c_str());
      if (first)
      {
        option = "PLC PMC";
      }
      else
      {
        option = "SAME PLC PMC";
      }
      first = false;
      th1DBG->DrawClone(option.c_str());
    }
  }
  else if (cname.find("VsPt") < cname.length())
  {
    for (auto etamin : list)
    {
      auto etamax = etamin + window;
      histo3D->GetYaxis()->SetRangeUser(etamin, etamax);
      std::string ytitle = "\\sigma (" + std::string(histo3D->GetZaxis()->GetTitle()) + ")";
      auto title = Form("_%1.2f_%1.2f_xz", etamin, etamax);
      auto aDBG = (TH2F *)histo3D->Project3D(title);
      aDBG->FitSlicesX(nullptr, 0, -1, 4, "QNR", &aSlices);
      auto th1DBG = (TH1F *)aSlices[iPar];
      th1DBG->SetTitle(Form("%1.2f < \\eta < %1.2f", etamin, etamax));
      th1DBG->SetStats(0);
      th1DBG->SetYTitle(ytitle.c_str());
      if (first)
      {
        option = "PLC PMC";
      }
      else
      {
        option = "SAME PLC PMC";
      }
      first = false;
      th1DBG->DrawClone(option.c_str());
    }
  }
  else
  {
    exit(1);
  }

  histo3D->GetYaxis()->SetRange(0, 0);
  histo3D->GetXaxis()->SetRange(0, 0);

  TPaveText *t = new TPaveText(0.2223748, 0.9069355, 0.7776252, 0.965, "brNDC"); // left-up
  t->SetBorderSize(0);
  t->SetFillColor(gStyle->GetTitleFillColor());
  t->AddText(ctitle.c_str());
  t->Draw();

  canvas->BuildLegend();
  canvas->SetTicky();
  canvas->SetGridy();
  if (0)
  {
    cname += ".png";
    canvas->Print(cname.c_str());
  }
}

//-------------------------------------------------------------------------------------------
void loadGlobalHistos()
{

  TObjArray *objar;

  TFile *f = new TFile(Form("GlobalForwardAssessment.root"));

  mTrackNumberOfClusters = std::unique_ptr<TH1F>((TH1F *)f->Get("mGlobalFwdNumberOfMFTClusters"));

  mTrackInvQPt = std::unique_ptr<TH1F>((TH1F *)f->Get("mGlobalFwdInvQPt"));

  mTrackChi2 = std::unique_ptr<TH1F>((TH1F *)f->Get("mGlobalFwdChi2"));

  mTrackCharge = std::unique_ptr<TH1F>((TH1F *)f->Get("mGlobalFwdCharge"));

  mTrackPhi = std::unique_ptr<TH1F>((TH1F *)f->Get("mGlobalFwdPhi"));

  mTrackEta = std::unique_ptr<TH1F>((TH1F *)f->Get("mGlobalFwdEta"));

  for (auto minNClusters : sMinNClustersList)
  {
    auto nHisto = minNClusters - sMinNClustersList[0];
    mTrackEtaNCls[nHisto] = std::unique_ptr<TH1F>((TH1F *)f->Get(Form("mGlobalFwdEta_%d_MinClusters", minNClusters)));

    mTrackPhiNCls[nHisto] = std::unique_ptr<TH1F>((TH1F *)f->Get(Form("mGlobalFwdPhi_%d_MinClusters", minNClusters)));

    mTrackXYNCls[nHisto] = std::unique_ptr<TH2F>((TH2F *)f->Get(Form("mGlobalFwdXY_%d_MinClusters", minNClusters)));

    mTrackEtaPhiNCls[nHisto] = std::unique_ptr<TH2F>((TH2F *)f->Get(Form("mGlobalFwdEtaPhi_%d_MinClusters", minNClusters)));
  }

  mTrackTanl = std::unique_ptr<TH1F>((TH1F *)f->Get("mGlobalFwdTanl"));

  // Creating MC-based histos
  if (mUseMC)
  {

    mHistPhiRecVsPhiGen = std::unique_ptr<TH2F>((TH2F *)f->Get("mGMTrackPhiRecVsPhiGen"));

    mHistEtaRecVsEtaGen = std::unique_ptr<TH2F>((TH2F *)f->Get("mGMTrackEtaRecVsEtaGen"));

    for (int trackType = 0; trackType < kNumberOfTrackTypes; trackType++)
    {
      mHistPhiVsEta[trackType] = std::unique_ptr<TH2F>((TH2F *)f->Get((std::string("mGMTrackPhiVsEta") + mNameOfTrackTypes[trackType]).c_str()));

      mHistPtVsEta[trackType] = std::unique_ptr<TH2F>((TH2F *)f->Get((std::string("mGMTrackPtVsEta") + mNameOfTrackTypes[trackType]).c_str()));

      mHistPhiVsPt[trackType] = std::unique_ptr<TH2F>((TH2F *)f->Get((std::string("mGMTrackPhiVsPt") + mNameOfTrackTypes[trackType]).c_str()));

      if (trackType != kReco)
      {
        mHistZvtxVsEta[trackType] = std::unique_ptr<TH2F>((TH2F *)f->Get((std::string("mGMTrackZvtxVsEta") + mNameOfTrackTypes[trackType]).c_str()));
      }
      if (trackType == kGen || trackType == kPairable)
      {
        mHistRVsZ[trackType] = std::unique_ptr<TH2F>((TH2F *)f->Get((std::string("mGMTrackRVsZ") + mNameOfTrackTypes[trackType]).c_str()));
      }
    }

    // Histos for Reconstruction assessment
    mChargeMatchEff = std::unique_ptr<TEfficiency>((TEfficiency *)f->Get("mGMTrackQMatchEff"));

    const int nTH3Histos = TH3Names.size();
    auto n3Histo = 0;
    for (auto &h : mTH3Histos)
    {
      h = std::unique_ptr<TH3F>((TH3F *)f->Get(TH3Names[n3Histo]));
      ++n3Histo;
    }
  }
}

//__________________________________________________________
void getHistos(TObjArray &objar)
{

  objar.Add(mTrackNumberOfClusters.get());
  objar.Add(mTrackInvQPt.get());
  objar.Add(mTrackChi2.get());
  objar.Add(mTrackCharge.get());
  objar.Add(mTrackPhi.get());
  objar.Add(mTrackEta.get());
  for (auto minNClusters : sMinNClustersList)
  {
    auto nHisto = minNClusters - sMinNClustersList[0];
    objar.Add(mTrackEtaNCls[nHisto].get());
    objar.Add(mTrackPhiNCls[nHisto].get());
    objar.Add(mTrackXYNCls[nHisto].get());
    objar.Add(mTrackEtaPhiNCls[nHisto].get());
  }
  objar.Add(mTrackTanl.get());

  if (mUseMC)
  {
    objar.Add(mHistPhiRecVsPhiGen.get());
    objar.Add(mHistEtaRecVsEtaGen.get());
    for (int TrackType = 0; TrackType < kNumberOfTrackTypes; TrackType++)
    {
      objar.Add(mHistPhiVsEta[TrackType].get());
      objar.Add(mHistPtVsEta[TrackType].get());
      objar.Add(mHistPhiVsPt[TrackType].get());
      objar.Add(mHistZvtxVsEta[TrackType].get());
      if (TrackType == kGen || TrackType == kPairable)
      {
        objar.Add(mHistRVsZ[TrackType].get());
      }
    }

    // Histos for Reconstruction assessment

    for (auto &h : mTH3Histos)
    {
      objar.Add(h.get());
    }

    objar.Add(mChargeMatchEff.get());
    objar.Add(mPurityPt.get());
    objar.Add(mPurityPtInner.get());
    objar.Add(mPurityPtOuter.get());
    objar.Add(mPairingEtaPt.get());

    if (1)
    {
      for (int slicedCanvas = 0; slicedCanvas < kNSlicedTH3; slicedCanvas++)
      {
        objar.Add(mSlicedCanvas[slicedCanvas]);
      }
      for (int matchingCanvas = 0; matchingCanvas < kNGMAssesmentCanvases; matchingCanvas++)
      {
        objar.Add(mAssessmentCanvas[matchingCanvas]);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------
void finalize()
{

  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0); // Remove title of first histogram from canvas
  gStyle->SetMarkerStyle(kFullCircle);
  gStyle->SetMarkerSize(1.5);

  std::vector<float> ptList({.5, 5., 10., 18.0});
  float ptWindow = 1.0;
  std::vector<float> etaList({2.5, 3.0});
  float etaWindow = 0.5;

  std::vector<float> sliceList;
  float sliceWindow;

  for (int nCanvas = 0; nCanvas < kNSlicedTH3; nCanvas++)
  {
    if (nCanvas % 2)
    {
      sliceList = etaList;
      sliceWindow = etaWindow;
    }
    else
    {
      sliceList = ptList;
      sliceWindow = ptWindow;
    }
    mSlicedCanvas[nCanvas] = new TCanvas(TH3SlicedNames[nCanvas], TH3SlicedNames[nCanvas], 1080, 1080);
    mSlicedCanvas[nCanvas]->UseCurrentStyle();
    mSlicedCanvas[nCanvas]->cd();
    TH3Slicer(mSlicedCanvas[nCanvas], mTH3Histos[TH3SlicedMap[nCanvas]], sliceList, sliceWindow, 2);
  }
  auto &Reco = mTH3Histos[kTH3GMTrackPtEtaMatchScore];
  auto &hTrue = mTH3Histos[kTH3GMTruePtEtaMatchScore];
  auto &hPairable = mTH3Histos[kTH3GMPairablePtEtaZ];

  auto RecoEtaPt = (TH2D *)Reco->Project3D("xy COLZ");
  auto PairableEtaPt = (TH2D *)hPairable->Project3D("xy COLZ");

  auto RecoPtProj = (TH1 *)Reco->ProjectionX();
  auto TruePtProj = (TH1 *)hTrue->ProjectionX();
  mPurityPt = std::make_unique<TEfficiency>(*TruePtProj, *RecoPtProj);
  mPurityPt->SetNameTitle("GMTrackPurity", "GMTrackPurity");
  auto minBin = Reco->GetYaxis()->FindBin(3.0);
  auto maxBin = Reco->GetYaxis()->FindBin(3.6);
  auto RecoPtProjInner = (TH1 *)Reco->ProjectionX("InnerReco", minBin, maxBin);
  auto TruePtProjInner = (TH1 *)hTrue->ProjectionX("InnerTrue", minBin, maxBin);
  mPurityPtInner = std::make_unique<TEfficiency>(*TruePtProjInner, *RecoPtProjInner);
  mPurityPtInner->SetNameTitle("GMTrackPurityInnerEta", "GMTrackPurity (3.0 < #eta < 3.6 )");

  minBin = Reco->GetYaxis()->FindBin(2.4);
  maxBin = Reco->GetYaxis()->FindBin(3.0);
  auto RecoPtProjOuter = (TH1 *)Reco->ProjectionX("OuterReco", minBin, maxBin);
  auto TruePtProjOuter = (TH1 *)hTrue->ProjectionX("OuterTrue", minBin, maxBin);
  mPurityPtOuter = std::make_unique<TEfficiency>(*TruePtProjOuter, *RecoPtProjOuter);
  mPurityPtOuter->SetNameTitle("GMTrackPurityOuterEta", "GMTrackPurity (2.4 < #eta < 3.0 )");

  mPairingEtaPt = (std::unique_ptr<TH2D>)static_cast<TH2D *>(RecoEtaPt->Clone());
  mPairingEtaPt->Divide(PairableEtaPt);
  mPairingEtaPt->SetNameTitle("GMTrackPairingEffEtaPt", "PairingEffEtaPt");
  mPairingEtaPt->SetOption("COLZ");
}

//-------------------------------------------------------------------------------------------
void finalize2()
{

  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0); // Remove title of first histogram from canvas
  gStyle->SetMarkerStyle(kFullCircle);
  gStyle->SetMarkerSize(1.0);

  auto &Reco = mTH3Histos[kTH3GMTrackPtEtaMatchScore];
  auto &hTrue = mTH3Histos[kTH3GMTruePtEtaMatchScore];
  auto &hPairable = mTH3Histos[kTH3GMPairablePtEtaZ];

  // Inner pseudorapidity
  auto minBin = Reco->GetYaxis()->FindBin(2.4);
  auto midBin = Reco->GetYaxis()->FindBin(3.0);
  auto maxBin = Reco->GetYaxis()->FindBin(3.6);
  auto PairablePtProjInner = (TH1 *)hTrue->ProjectionX("PairableInner", midBin, maxBin);
  auto PairablePtProjOuter = (TH1 *)hTrue->ProjectionX("PairableOuter", minBin, midBin);

  auto RecoEtaPt = (TH2D *)Reco->Project3D("xy COLZ");
  auto PairableEtaPt = (TH2D *)hPairable->Project3D("xy COLZ");
  auto PairablePt = (TH1D *)hPairable->Project3D("x");

  /// Purity vs score cuts
  float maxCut = 15.f;
  int nSteps = 15;
  float scoreStep = maxCut / nSteps;
  for (float scoreCut = scoreStep; scoreCut < maxCut; scoreCut += scoreStep)
  {
    std::cout << "scoreCut = " << scoreCut << " ===> " << Form("%.2f", scoreCut) << std::endl;

    auto RecoPtProj = (TH1 *)Reco->ProjectionX(Form("_RecoPtProj%.2f", scoreCut));
    auto TruePtProj = (TH1 *)hTrue->ProjectionX(Form("_TruePtProj%.2f", scoreCut));
    mPurityPt = std::make_unique<TEfficiency>(*TruePtProj, *RecoPtProj);
    mPurityPt->SetNameTitle(Form("GMTrackPurityCut_%.2f", scoreCut), Form("GMTrackPurityCut_%.2f", scoreCut));

    // Inner pseudorapidity
    auto maxScoreBin = Reco->GetZaxis()->FindBin(scoreCut);
    auto RecoPtProjInner = (TH1 *)Reco->ProjectionX(Form("_InnerRecoCut_%.2f", scoreCut), midBin, maxBin, 0, maxScoreBin);
    auto TruePtProjInner = (TH1 *)hTrue->ProjectionX(Form("_InnerTrueCut_%.2f", scoreCut), midBin, maxBin, 0, maxScoreBin);

    mPurityPtInnerVec.emplace_back(std::make_unique<TEfficiency>(*TruePtProjInner, *RecoPtProjInner));
    mPurityPtInnerVec.back()->SetNameTitle(Form("GMTrackPurityInnerEtaCut_%.2f", scoreCut), Form("GMTrackPurity (3.0 < #eta < 3.6 ) cut %.2f", scoreCut));

    auto &hPInner = mPurityPtInnerVecTH2.emplace_back((std::unique_ptr<TH2D>)static_cast<TH2D *>(TruePtProjInner->Clone()));
    hPInner->Divide(RecoPtProjInner);
    hPInner->SetNameTitle(Form("TH2GMTrackPurityInnerEtaCut_%.2f", scoreCut), Form("%.2f cut", scoreCut));
    hPInner->SetOption("COLZ");
    hPInner->SetMarkerStyle(kFullCircle);
    hPInner->SetMinimum(0.0);
    hPInner->SetMaximum(1.2);

    auto &hInner = mPairingPtInnerVecTH1.emplace_back((std::unique_ptr<TH1D>)static_cast<TH1D *>(RecoPtProjInner->Clone()));
    hInner->Divide(PairablePtProjInner);
    hInner->SetNameTitle(Form("GMTrackPairingEffInnerPtCut_%.2f", scoreCut), Form("%.2f cut", scoreCut));
    hInner->SetOption("COLZ");
    hInner->SetMarkerStyle(kFullCircle);
    hInner->SetMinimum(0.0);
    hInner->SetMaximum(1.8);

    // Outer pseudorapidity
    auto RecoPtProjOuter = (TH1 *)Reco->ProjectionX(Form("_OuterRecoCut_%.2f", scoreCut), minBin, midBin, 0, maxScoreBin);
    auto TruePtProjOuter = (TH1 *)hTrue->ProjectionX(Form("_OuterTrueCut_%.2f", scoreCut), minBin, midBin, 0, maxScoreBin);

    mPurityPtOuterVec.emplace_back(std::make_unique<TEfficiency>(*TruePtProjOuter, *RecoPtProjOuter));
    mPurityPtOuterVec.back()->SetNameTitle(Form("GMTrackPurityOuterEtaCut_%.2f", scoreCut), Form("GMTrackPurity (2.4 < #eta < 3.0 ) cut %.2f", scoreCut));

    auto &hPOuter = mPurityPtOuterVecTH2.emplace_back((std::unique_ptr<TH2D>)static_cast<TH2D *>(TruePtProjOuter->Clone()));
    hPOuter->Divide(RecoPtProjOuter);
    hPOuter->SetNameTitle(Form("TH2GMTrackPurityOuterEtaCut_%.2f", scoreCut), Form("%.2f cut", scoreCut));
    hPOuter->SetOption("COLZ");
    hPOuter->SetMarkerStyle(kFullCircle);
    hPOuter->SetMinimum(0.0);
    hPOuter->SetMaximum(1.2);

    auto &hOuter = mPairingPtOuterVecTH1.emplace_back((std::unique_ptr<TH1D>)static_cast<TH1D *>(RecoPtProjOuter->Clone()));
    hOuter->Divide(PairablePtProjInner);
    hOuter->SetNameTitle(Form("GMTrackPairingEffOuterPtCut_%.2f", scoreCut), Form("%.2f cut", scoreCut));
    hOuter->SetOption("COLZ");
    hOuter->SetMarkerStyle(kFullCircle);
    hOuter->SetMinimum(0.0);
    hOuter->SetMaximum(1.8);

    mPairingEtaPtVec.emplace_back((std::unique_ptr<TH2D>)static_cast<TH2D *>(RecoEtaPt->Clone()));
    mPairingEtaPtVec.back()->Divide(PairableEtaPt);
    mPairingEtaPtVec.back()->SetNameTitle(Form("GMTrackPairingEffEtaPtCut_%.2f", scoreCut), Form("%.2f", scoreCut));
    mPairingEtaPtVec.back()->SetOption("COLZ");
  }

  auto nCanvas = kPurityPtOuter;
  auto canvasName = GMAssesmentNames[nCanvas];
  mAssessmentCanvas[nCanvas] = new TCanvas(canvasName, canvasName, 1080, 1080);
  mAssessmentCanvas[nCanvas]->UseCurrentStyle();
  mAssessmentCanvas[nCanvas]->cd();
  //auto canvas = new TCanvas("PurityPtOuter", "PurityPtOuter", 1024, 800);
  auto first = true;
  std::string option;

  for (auto &th2 : mPurityPtOuterVecTH2)
  {
    if (first)
    {
      option = "hist P PMC";
    }
    else
    {
      option = "hist SAME P PMC";
    }
    first = false;
    th2->Draw(option.c_str());
  }
  TPaveText *t = new TPaveText(0.2223748, 0.9069355, 0.7776252, 0.965, "brNDC"); // left-up
  t->SetBorderSize(0);
  t->SetFillColor(gStyle->GetTitleFillColor());
  t->AddText("Global Muon Track Purity (2.4 < #eta < 3.0 )");
  t->Draw();

  mAssessmentCanvas[nCanvas]->BuildLegend(.8, .15, .96, .87);
  mAssessmentCanvas[nCanvas]->SetTicky();
  mAssessmentCanvas[nCanvas]->SetGridy();

  nCanvas = kPurityPtInner;
  canvasName = GMAssesmentNames[nCanvas];
  mAssessmentCanvas[nCanvas] = new TCanvas(canvasName, canvasName, 1080, 1080);
  mAssessmentCanvas[nCanvas]->UseCurrentStyle();
  mAssessmentCanvas[nCanvas]->cd();
  first = true;

  for (auto &th2 : mPurityPtInnerVecTH2)
  {
    if (first)
    {
      option = "hist P PMC";
    }
    else
    {
      option = "hist SAME P PMC";
    }
    first = false;
    th2->Draw(option.c_str());
  }
  t = new TPaveText(0.2223748, 0.9069355, 0.7776252, 0.965, "brNDC");
  t->SetBorderSize(0);
  t->SetFillColor(gStyle->GetTitleFillColor());
  t->AddText("Global Muon Track Purity (3.0 < #eta < 3.6 )");
  t->Draw();

  mAssessmentCanvas[nCanvas]->BuildLegend(.8, .15, .96, .87);
  mAssessmentCanvas[nCanvas]->SetTicky();
  mAssessmentCanvas[nCanvas]->SetGridy();

  nCanvas = kPairingEffPtOuter;
  canvasName = GMAssesmentNames[nCanvas];
  mAssessmentCanvas[nCanvas] = new TCanvas(canvasName, canvasName, 1080, 1080);
  mAssessmentCanvas[nCanvas]->UseCurrentStyle();
  mAssessmentCanvas[nCanvas]->cd();
  first = true;

  for (auto &th2 : mPairingPtOuterVecTH1)
  {
    if (first)
    {
      option = "hist P PMC";
    }
    else
    {
      option = "hist SAME P PMC";
    }
    first = false;
    th2->Draw(option.c_str());
  }
  t = new TPaveText(0.2223748, 0.9069355, 0.7776252, 0.965, "brNDC");
  t->SetBorderSize(0);
  t->SetFillColor(gStyle->GetTitleFillColor());
  t->AddText("Global Muon Track Pairing Efficiency (3.0 < #eta < 3.6 )");
  t->Draw();

  mAssessmentCanvas[nCanvas]->BuildLegend(.8, .15, .96, .87);
  mAssessmentCanvas[nCanvas]->SetTicky();
  mAssessmentCanvas[nCanvas]->SetGridy();

  nCanvas = kPairingEffPtInner;
  canvasName = GMAssesmentNames[nCanvas];
  mAssessmentCanvas[nCanvas] = new TCanvas(canvasName, canvasName, 1080, 1080);
  mAssessmentCanvas[nCanvas]->UseCurrentStyle();
  mAssessmentCanvas[nCanvas]->cd();
  first = true;

  for (auto &th2 : mPairingPtInnerVecTH1)
  {
    if (first)
    {
      option = "hist P PMC";
    }
    else
    {
      option = "hist SAME P PMC";
    }
    first = false;
    th2->Draw(option.c_str());
  }
  t = new TPaveText(0.2223748, 0.9069355, 0.7776252, 0.965, "brNDC");
  t->SetBorderSize(0);
  t->SetFillColor(gStyle->GetTitleFillColor());
  t->AddText("Global Muon Track Pairing Efficiency (2.4 < #eta < 3.0 )");
  t->Draw();

  mAssessmentCanvas[nCanvas]->BuildLegend(.8, .15, .96, .87);
  mAssessmentCanvas[nCanvas]->SetTicky();
  mAssessmentCanvas[nCanvas]->SetGridy();
}

//-------------------------------------------------------------------------------------------
void FwdAssessmentFinalization()
{
  loadGlobalHistos();
  finalize2();

  TFile *fout = new TFile("GlobalFwdAssessmentFinalized.root", "RECREATE");
  TObjArray objarOut;
  getHistos(objarOut);
  objarOut.Write();
  fout->Close();
}