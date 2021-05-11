
#ifndef MUON_MATCHING
#define MUON_MATCHING

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsPassive/Absorber.h"
#include "DetectorsPassive/Cave.h"
#include "DetectorsPassive/Shil.h"
#include "Field/MagneticField.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTTracking/IOUtils.h"
#include "MFTTracking/TrackCA.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TGeoManager.h"
#include <TGeoGlobalMagField.h>

#include "Math/SMatrix.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"
#include "TStopwatch.h"
#endif

#include "MCHTracking/TrackParam.h"
#include "MFTBase/Constants.h"
#include "MFTTracking/Cluster.h"
#include "MFTTracking/Constants.h"
#include "include/GlobalMuonTrack.h"
#include "include/TrackExtrap.h"
#include "include/tempMCHTrack.h"
#include "include/MLHelpers.h"

using namespace TMVA;

#include <iostream>

using MCHTrack = o2::mch::TrackParam;
using MFTTrack = o2::mft::TrackMFT;
using MCHTrackConv = o2::track::GlobalMuonTrack; // MCHTracks converted to the Fwd coordinate system

using GlobalMuonTrack = o2::track::GlobalMuonTrack;
using GlobalMuonTrackExt = o2::track::GlobalMuonTrackExt;
using MCLabels = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
using MFTCluster = o2::mft::Cluster;

using SMatrix22 = ROOT::Math::SMatrix<double, 2>;
using SMatrix25 = ROOT::Math::SMatrix<double, 2, 5>;
using SMatrix52 = ROOT::Math::SMatrix<double, 5, 2>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using SVector2 = ROOT::Math::SVector<double, 2>;
using SVector4 = ROOT::Math::SVector<double, 4>;
using SVector5 = ROOT::Math::SVector<double, 5>;

using SMatrix44 = ROOT::Math::SMatrix<double, 4>;
using SMatrix45 = ROOT::Math::SMatrix<double, 4, 5>;
using SMatrix54 = ROOT::Math::SMatrix<double, 5, 4>;
using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
using SMatrix55Sym =
  ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

//_________________________________________________________________________________________________
class MUONMatcher
{
 public:
  MUONMatcher();
  ~MUONMatcher() = default;
  void Clear();
  void SetMatchingPlaneZ(double z) {
    mMatchingPlaneZ = z;

    double SAbsZBeg =  -90.; ///< Position of the begining of the absorber (cm)
    double SAbsZEnd = -505.; ///< Position of the end of the absorber (cm)

    // Check the matching plane position with respect to the absorber (spectro
    // z<0)
    if (mMatchingPlaneZ < SAbsZBeg) {
      if (mMatchingPlaneZ < SAbsZEnd) {
	LOG(WARNING) << "Ending Z (" << mMatchingPlaneZ
		     << ") downstream the front absorber (zAbsorberEnd = "
		     << SAbsZEnd << ")";
      }
      else {
	LOG(WARNING) << "Ending Z (" << mMatchingPlaneZ
		     << ") inside the front absorber (" << SAbsZBeg << ", "
		     << SAbsZEnd << ")";
      }
    }
  }
  void SetVerbosity(bool v = true) { mVerbose = v; }
  void LoadAbsorber();

  // Track IO
  void loadMCHTracks();
  void loadMFTTracksOut();
  void saveGlobalMuonTracks();

  void printMFTLabels()
  {
    auto i = 0;
    for (auto label : mftTrackLabels) {
      std::cout << " Track " << i << " label: ";
      label.print();
      std::cout << std::endl;
      i++;
    }
  }

  void initGlobalTracks(); // Configure Global Tracks with MCH track parameters
  void fitTracks();          // Fit all matched tracks

  // Matching methods
  // Position
  double
    matchMFT_MCH_TracksXY(const MCHTrackConv& mchTrack,
                          const MFTTrack& mftTrack); // Compute track matching
  //// Position & Angles
  double matchMFT_MCH_TracksXYPhiTanl(const MCHTrackConv& mchTrack,
                                      const MFTTrack& mftTrack);
  //// Position, Angles & Charged Momentum
  double matchMFT_MCH_TracksAllParam(const MCHTrackConv& mchTrack,
                                     const MFTTrack& mftTrack);
  //// Hiroshima's Matching
  double matchHiroshima(const GlobalMuonTrack& mchTrack,
                        const MFTTrack& mftTrack);
  //// Matching using trained ML
  double matchTrainedML(const MCHTrackConv& mchTrack,
                        const MFTTrack& mftTrack);

  void setMatchingFunction(double (MUONMatcher::*func)(const MCHTrackConv&,
                                                       const MFTTrack&))
  {
    mMatchFunc = func;
    if (func == &MUONMatcher::matchMFT_MCH_TracksXY) {
      mMatchingHelper.MatchingFunction = "_matchXY";
    }
    if (func == &MUONMatcher::matchMFT_MCH_TracksXYPhiTanl) {
      mMatchingHelper.MatchingFunction = "_matchXYPhiTanl";
    }
    if (func == &MUONMatcher::matchMFT_MCH_TracksAllParam) {
      mMatchingHelper.MatchingFunction = "_matchAllParams";
    }
    if (func == &MUONMatcher::matchHiroshima) {
      mMatchingHelper.MatchingFunction = "_matchHiroshima";
    }
    if (func == &MUONMatcher::matchTrainedML) {
      mMatchingHelper.MatchingFunction = "_matchML";
    } else { // Clear MLFeaturesFunction string if not using ML
      mMatchingHelper.MLFeaturesFunction = "";
    }

    std::cout << " ** MUONMATCHER: Setting matching function => "
              << mMatchingHelper.MatchingFunction << std::endl;
  }
  void setCustomMatchingFunction(double (*func)(const MCHTrackConv&,
                                                const MFTTrack&),
                                 std::string nickname)
  {
    mCustomMatchFunc = func;
    mMatchingHelper.MatchingFunction = nickname;
    std::cout << " ** MUONMATCHER: Setting custom matching function => "
              << mMatchingHelper.MatchingFunction << std::endl;
  }
  void setMLFeatureFunction(void (*func)(const MCHTrackConv&,
                                         const MFTTrack&,
                                         float*),
                            int nFeatures, std::string nickname,
                            void (*funcNames)(string*) = nullptr)
  {
    mMLFeaturesFunc = func;
    mNInputFeatures = nFeatures;
    mMatchingHelper.MLFeaturesFunction = nickname;
    if (funcNames) {
      (*funcNames)(mMLInputFeaturesName);
    } else {
      std::cout << "No ML input feature names informed. Setting defaults Feature_N." << std::endl;
      for (std::size_t nFeature = 0; nFeature < mNInputFeatures; nFeature++) {
        mMLInputFeaturesName[nFeature] = Form("Feature_%d", (Int_t)nFeature);
      }
    }
  }
  void setMatchSaveAll(bool val) { mMatchSaveAll = val; }
  void
    runHeavyMatching();    // Finds best match (no search cut, no event separation)
  void runEventMatching(); // Finds best match event-per-event
  bool printMatchingPlaneView(int event, int MCHTrackID = 0);
  void exportNMatchingPlaneViews(int nTracks)
  {
    if (mMatchSaveAll) {
      std::cout << " ** exportNMatchingPlaneViews disabled for option --matchSaveAll" << std::endl;
      return;
    }
    loadMCHTracks();
    initGlobalTracks();
    auto event = 0;
    if (nTracks > mMatchingHelper.nMCHTracks)
      nTracks = mMatchingHelper.nMCHTracks;
    std::cout << " ** exportNMatchingPlaneViews: " << nTracks << std::endl;
    for (auto& mchTracks : mSortedGlobalMuonTracks) {
      auto MCHTrackID = 0;
      for (auto& mchTrack : mchTracks) {
        if (nTracks > 0)
          printMatchingPlaneView(event, MCHTrackID);
        nTracks--;
        MCHTrackID++;
      }
      event++;
    }
  };

  void exportTrainingDataRoot(int nMCHTracks = -1);
  void setCorrectMatchIgnoreCut(bool v = true) { mCorrectMatchIgnoreCut = v; };

  // Matching cuts
  void enableChargeMatchCut() { mChargeCutEnabled = true; }
  bool matchingCut(const GlobalMuonTrack&,
                   const MFTTrack&); // Calls configured cut function

  double matchingEval(const MCHTrackConv&,
                      const MFTTrack&); // Calls configured cut function
  void setCutFunction(bool (MUONMatcher::*func)(const MCHTrackConv&,
                                                const MFTTrack&));
  void setCustomCutFunction(bool (*func)(const MCHTrackConv&,
                                         const MFTTrack&))
  {
    mCustomCutFunc = func;
  }

  // Machine Learning Methods
  void EvaluateML();
  void setMLFeatures(const MCHTrackConv& mchTr, const MFTTrack& mftTr) { (*mMLFeaturesFunc)(mchTr, mftTr, &mMLInputFeatures[0]); }

  void configureTMVA(std::string filename, float scorecut)
  {
    mTMVAWeightFileName = filename;
    mMLScoreCut = scorecut;
    mTMVAReader = new TMVA::Reader("!Color:!Silent");
    for (std::size_t i = 0; i < mNInputFeatures; i++) {
      mTMVAReader->AddVariable(mMLInputFeaturesName[i], &mMLInputFeatures[i]);
    }

    mTMVAReader->BookMVA("MUONMatcherML", mTMVAWeightFileName);
  }
  void MLRegression(std::string input_name, std::string trainingfile,
                    std::string trainingstr);
  void MLClassification(std::string input_name, std::string trainingfile,
                        std::string trainingstr);

  void MLTraining()
  {
    std::string network_ID;
    std::string MLLayout("");
    std::string MLStrat("");
    std::string MLOpt("");
    std::string training_file = gSystem->Getenv("ML_TRAINING_FILE");

    if (gSystem->Getenv("ML_LAYOUT")) {
      MLLayout = gSystem->Getenv("ML_LAYOUT");
      network_ID += MLLayout + "_";
    }
    if (gSystem->Getenv("ML_TRAINING_STRAT")) {
      MLStrat = gSystem->Getenv(
        "ML_TRAINING_STRAT"); //#TODO reorganize these strings in MLHelpers
      network_ID += MLStrat + "_";
    }
    if (gSystem->Getenv("ML_GENERAL_OPT")) {
      MLOpt = gSystem->Getenv("ML_GENERAL_OPT");
      network_ID += MLOpt + "_";
    }

    std::string training_string("");
    if (network_ID != "") {
      training_string = opt_reader();
    } else {
      std::cout << " [WARNING] Configurations for ML method were not setted. I'll use TMVAs default; hope that works!" << endl;
    }

    std::cout << " Network name: " << network_ID << "\n"
              << std::endl;
    std::cout << " Complete Option String: " << training_string << "\n"
              << std::endl;

    std::string MLAnalysisType = gSystem->Getenv("ML_TYPE");
    if (MLAnalysisType == "regression" || MLAnalysisType == "Regression") {
      MLRegression(network_ID, training_file, training_string);

    } else if (MLAnalysisType == "Classification" || MLAnalysisType == "classification") {
      MLClassification(network_ID, training_file, training_string);

    } else {
      std::cout << " Type of ML analysis does not exist or it is not supported right now. Please choose Classification or Regression. " << std::endl;
    }
  }

  //  Built-in cut functions
  bool matchCutDisabled(const MCHTrackConv&, const MFTTrack&);
  bool matchCutDistance(const MCHTrackConv&, const MFTTrack&);
  bool matchCutDistanceAndAngles(const MCHTrackConv&, const MFTTrack&);
  bool matchCutDistanceSigma(const MCHTrackConv&, const MFTTrack&);
  bool matchCut3SigmaXYAngles(const MCHTrackConv&, const MFTTrack&);
  bool matchCutVarXYAngles(const MCHTrackConv&, const MFTTrack&);
  void setCutParam(int index, double param)
  {
    if (index > ((int)mCutParams.size() - 1))
      mCutParams.resize(index + 1);
    mCutParams[index] = param;
    std::cout << " ** MUONMATCHER: Setting matching cutParam[" << index
              << "] = " << param << std::endl;
  };

 private:
  // Private IO methods
  void loadMFTClusters();

  // Track methods
  GlobalMuonTrack
    MCHtoGlobal(MCHTrack&); // Convert MCH Track to GlobalMuonTrack;
  void finalize();

  // Global Muon Track Methods
  void fitGlobalMuonTrack(GlobalMuonTrack&); // Kalman filter
  bool computeCluster(GlobalMuonTrack&, MFTCluster&);
  double (MUONMatcher::*mMatchFunc)(const MCHTrackConv&, const MFTTrack&);
  double (*mCustomMatchFunc)(const MCHTrackConv&,
                             const MFTTrack&) = nullptr;
  bool (MUONMatcher::*mCutFunc)(const MCHTrackConv&, const MFTTrack&);
  bool (*mCustomCutFunc)(const GlobalMuonTrack&, const MFTTrack&) = nullptr;

  void (*mMLFeaturesFunc)(const MCHTrackConv&, const MFTTrack&, float*);

  // Data Members
  std::vector<MFTTrack> mMFTTracks;
  std::vector<std::vector<MFTTrack>> mSortedMFTTracks;

  std::vector<std::list<MCHTrack>> mSortedMCHTracks;

  std::vector<std::vector<GlobalMuonTrack>> mSortedGlobalMuonTracks;

  std::vector<std::vector<GlobalMuonTrackExt>> mSortedGlobalMuonTracksExt;

  std::vector<MFTCluster> mMFTClusters;
  std::vector<int> mtrackExtClsIDs;
  std::vector<o2::itsmft::ROFRecord> mMFTTracksROFs;
  int mNEvents = 0;
  MatchingHelper mMatchingHelper;

  std::vector<o2::MCCompLabel> mftTrackLabels;
  std::vector<std::vector<int>> mftTrackLabelsIDx;

  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>> mSortedMCHTrackLabels;

  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>> mSortedGlobalTrackLabels;

  // MCH Track Propagation clasee
  o2::mch::TrackExtrap mMCHTrackExtrap;

  double mField_z;
  const double sLastMFTPlaneZ = -77.5;
  double mMatchingPlaneZ = sLastMFTPlaneZ;
  std::vector<double> mCutParams;
  bool mVerbose = false;
  TGeoManager* mGeoManager;
  bool mMatchSaveAll = false;
  bool mChargeCutEnabled = false;

  // TMVA interface
  static const std::size_t sMaxMLFeatures = 50;
  std::string mTMVAWeightFileName;
  TMVA::Reader* mTMVAReader;
  float_t mMLScoreCut;
  float_t mMLInputFeatures[sMaxMLFeatures];
  std::size_t mNInputFeatures = -1;
  bool mCorrectMatchIgnoreCut = false; // Cut function not applied to correct match on training data

  string mMLInputFeaturesName[sMaxMLFeatures];
};

//_________________________________________________________________________________________________
std::string getParamString(o2::track::TrackParCovFwd t)
{
  std::string param;
  param = " x = " + std::to_string(t.getX()) +
          " y = " + std::to_string(t.getY()) +
          " phi = " + std::to_string(t.getPhi()) +
          " tanl = " + std::to_string(t.getTanl()) +
          " q*pt = " + std::to_string(1.0 / t.getInvQPt());
  return param;
}

//_________________________________________________________________________________________________
std::string getCovString(o2::track::TrackParCovFwd t)
{
  std::string param;
  param = "Cov: (" + std::to_string(t.getCovariances()(0, 0)) + " ; " +
          std::to_string(t.getCovariances()(1, 1)) + " ; " +
          std::to_string(t.getCovariances()(2, 2)) + " ; " +
          std::to_string(t.getCovariances()(3, 3)) + " ; " +
          std::to_string(t.getCovariances()(4, 4)) + ")";
  return param;
}

//_________________________________________________________________________________________________
void MLParCov40Features(const MCHTrackConv& mchTrack, const MFTTrack& mftTrack, float* features)
{

  features[0] = mftTrack.getX();
  features[1] = mftTrack.getY();
  features[2] = mftTrack.getPhi(),
  features[3] = mftTrack.getTanl();
  features[4] = mftTrack.getInvQPt();
  features[5] = mftTrack.getCovariances()(0, 0);
  features[6] = mftTrack.getCovariances()(0, 1);
  features[7] = mftTrack.getCovariances()(1, 1);
  features[8] = mftTrack.getCovariances()(0, 2);
  features[9] = mftTrack.getCovariances()(1, 2);
  features[10] = mftTrack.getCovariances()(2, 2);
  features[11] = mftTrack.getCovariances()(0, 3);
  features[12] = mftTrack.getCovariances()(1, 3);
  features[13] = mftTrack.getCovariances()(2, 3);
  features[14] = mftTrack.getCovariances()(3, 3);
  features[15] = mftTrack.getCovariances()(0, 4);
  features[16] = mftTrack.getCovariances()(1, 4);
  features[17] = mftTrack.getCovariances()(2, 4);
  features[18] = mftTrack.getCovariances()(3, 4);
  features[19] = mftTrack.getCovariances()(4, 4);
  features[20] = mchTrack.getX();
  features[21] = mchTrack.getY();
  features[22] = mchTrack.getPhi(),
  features[23] = mchTrack.getTanl();
  features[24] = mchTrack.getInvQPt();
  features[25] = mchTrack.getCovariances()(0, 0);
  features[26] = mchTrack.getCovariances()(0, 1);
  features[27] = mchTrack.getCovariances()(1, 1);
  features[28] = mchTrack.getCovariances()(0, 2);
  features[29] = mchTrack.getCovariances()(1, 2);
  features[30] = mchTrack.getCovariances()(2, 2);
  features[31] = mchTrack.getCovariances()(0, 3);
  features[32] = mchTrack.getCovariances()(1, 3);
  features[33] = mchTrack.getCovariances()(2, 3);
  features[34] = mchTrack.getCovariances()(3, 3);
  features[35] = mchTrack.getCovariances()(0, 4);
  features[36] = mchTrack.getCovariances()(1, 4);
  features[37] = mchTrack.getCovariances()(2, 4);
  features[38] = mchTrack.getCovariances()(3, 4);
  features[39] = mchTrack.getCovariances()(4, 4);
}

//_________________________________________________________________________________________________
void MLParCov40FeaturesNames(string* featuresNames)
{
  std::cout << "Setting features names: MLParCov40FeaturesNames" << std::endl;
  featuresNames[0] = "MFT_X";
  featuresNames[1] = "MFT_Y";
  featuresNames[2] = "MFT_Phi";
  featuresNames[3] = "MFT_Tanl";
  featuresNames[4] = "MFT_InvQPt";
  featuresNames[5] = "MFT_Cov00";
  featuresNames[6] = "MFT_Cov01";
  featuresNames[7] = "MFT_Cov11";
  featuresNames[8] = "MFT_Cov02";
  featuresNames[9] = "MFT_Cov12";
  featuresNames[10] = "MFT_Cov22";
  featuresNames[11] = "MFT_Cov03";
  featuresNames[12] = "MFT_Cov13";
  featuresNames[13] = "MFT_Cov23";
  featuresNames[14] = "MFT_Cov33";
  featuresNames[15] = "MFT_Cov04";
  featuresNames[16] = "MFT_Cov14";
  featuresNames[17] = "MFT_Cov24";
  featuresNames[18] = "MFT_Cov34";
  featuresNames[19] = "MFT_Cov44";
  featuresNames[20] = "MCH_X";
  featuresNames[21] = "MCH_Y";
  featuresNames[22] = "MCH_Phi";
  featuresNames[23] = "MCH_Tanl";
  featuresNames[24] = "MCH_InvQPt";
  featuresNames[25] = "MCH_Cov00";
  featuresNames[26] = "MCH_Cov01";
  featuresNames[27] = "MCH_Cov11";
  featuresNames[28] = "MCH_Cov02";
  featuresNames[29] = "MCH_Cov12";
  featuresNames[30] = "MCH_Cov22";
  featuresNames[31] = "MCH_Cov03";
  featuresNames[32] = "MCH_Cov13";
  featuresNames[33] = "MCH_Cov23";
  featuresNames[34] = "MCH_Cov33";
  featuresNames[35] = "MCH_Cov04";
  featuresNames[36] = "MCH_Cov14";
  featuresNames[37] = "MCH_Cov24";
  featuresNames[38] = "MCH_Cov34";
  featuresNames[39] = "MCH_Cov44";
}

//_________________________________________________________________________________________________
void MLParCovChiNPts42Features(const MCHTrackConv& mchTrack, const MFTTrack& mftTrack, float* features)
{

  features[0] = mftTrack.getX();
  features[1] = mftTrack.getY();
  features[2] = mftTrack.getPhi(),
  features[3] = mftTrack.getTanl();
  features[4] = mftTrack.getInvQPt();
  features[5] = mftTrack.getCovariances()(0, 0);
  features[6] = mftTrack.getCovariances()(0, 1);
  features[7] = mftTrack.getCovariances()(1, 1);
  features[8] = mftTrack.getCovariances()(0, 2);
  features[9] = mftTrack.getCovariances()(1, 2);
  features[10] = mftTrack.getCovariances()(2, 2);
  features[11] = mftTrack.getCovariances()(0, 3);
  features[12] = mftTrack.getCovariances()(1, 3);
  features[13] = mftTrack.getCovariances()(2, 3);
  features[14] = mftTrack.getCovariances()(3, 3);
  features[15] = mftTrack.getCovariances()(0, 4);
  features[16] = mftTrack.getCovariances()(1, 4);
  features[17] = mftTrack.getCovariances()(2, 4);
  features[18] = mftTrack.getCovariances()(3, 4);
  features[19] = mftTrack.getCovariances()(4, 4);

  features[20] = mchTrack.getX();
  features[21] = mchTrack.getY();
  features[22] = mchTrack.getPhi(),
  features[23] = mchTrack.getTanl();
  features[24] = mchTrack.getInvQPt();
  features[25] = mchTrack.getCovariances()(0, 0);
  features[26] = mchTrack.getCovariances()(0, 1);
  features[27] = mchTrack.getCovariances()(1, 1);
  features[28] = mchTrack.getCovariances()(0, 2);
  features[29] = mchTrack.getCovariances()(1, 2);
  features[30] = mchTrack.getCovariances()(2, 2);
  features[31] = mchTrack.getCovariances()(0, 3);
  features[32] = mchTrack.getCovariances()(1, 3);
  features[33] = mchTrack.getCovariances()(2, 3);
  features[34] = mchTrack.getCovariances()(3, 3);
  features[35] = mchTrack.getCovariances()(0, 4);
  features[36] = mchTrack.getCovariances()(1, 4);
  features[37] = mchTrack.getCovariances()(2, 4);
  features[38] = mchTrack.getCovariances()(3, 4);
  features[39] = mchTrack.getCovariances()(4, 4);

  features[40] = mftTrack.getTrackChi2();
  features[41] = mftTrack.getNumberOfPoints();
}

//_________________________________________________________________________________________________
void MLParCovChiNPts42FeaturesNames(string* featuresNames)
{
  std::cout << "Setting features names: MLParCovChiNPts42FeaturesNames" << std::endl;

  featuresNames[0] = "MFT_X";
  featuresNames[1] = "MFT_Y";
  featuresNames[2] = "MFT_Phi";
  featuresNames[3] = "MFT_Tanl";
  featuresNames[4] = "MFT_InvQPt";
  featuresNames[5] = "MFT_Cov00";
  featuresNames[6] = "MFT_Cov01";
  featuresNames[7] = "MFT_Cov11";
  featuresNames[8] = "MFT_Cov02";
  featuresNames[9] = "MFT_Cov12";
  featuresNames[10] = "MFT_Cov22";
  featuresNames[11] = "MFT_Cov03";
  featuresNames[12] = "MFT_Cov13";
  featuresNames[13] = "MFT_Cov23";
  featuresNames[14] = "MFT_Cov33";
  featuresNames[15] = "MFT_Cov04";
  featuresNames[16] = "MFT_Cov14";
  featuresNames[17] = "MFT_Cov24";
  featuresNames[18] = "MFT_Cov34";
  featuresNames[19] = "MFT_Cov44";

  featuresNames[20] = "MCH_X";
  featuresNames[21] = "MCH_Y";
  featuresNames[22] = "MCH_Phi";
  featuresNames[23] = "MCH_Tanl";
  featuresNames[24] = "MCH_InvQPt";
  featuresNames[25] = "MCH_Cov00";
  featuresNames[26] = "MCH_Cov01";
  featuresNames[27] = "MCH_Cov11";
  featuresNames[28] = "MCH_Cov02";
  featuresNames[29] = "MCH_Cov12";
  featuresNames[30] = "MCH_Cov22";
  featuresNames[31] = "MCH_Cov03";
  featuresNames[32] = "MCH_Cov13";
  featuresNames[33] = "MCH_Cov23";
  featuresNames[34] = "MCH_Cov33";
  featuresNames[35] = "MCH_Cov04";
  featuresNames[36] = "MCH_Cov14";
  featuresNames[37] = "MCH_Cov24";
  featuresNames[38] = "MCH_Cov34";
  featuresNames[39] = "MCH_Cov44";

  featuresNames[40] = "MFT_TrackChi2";
  featuresNames[41] = "MFT_NClust";
}

#include "MUONMatcher.cxx"

#endif /* MUON_MATCHING */
