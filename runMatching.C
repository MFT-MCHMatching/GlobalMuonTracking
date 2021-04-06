#if !defined(__CLING__) || defined(__ROOTCLING__)

#endif

#include "MUONMatcher.h"

//#ifdef __MAKECINT__
#pragma link C++ class GlobalMuonTrack + ;
#pragma link C++ class GlobalMuonTrackExt + ;
#pragma link C++ class std::vector < GlobalMuonTrack> + ;
#pragma link C++ class std::vector < GlobalMuonTrackExt> + ;
#pragma link C++ class tempMCHTrack + ;
#pragma link C++ class o2::mch::TrackExtrap + ;
#pragma link C++ class std::vector < tempMCHTrack> + ;
//#endif

MUONMatcher matcher;

//_________________________________________________________________________________________________
// Sample custom matching function that can be passed to MUONMatcher
double MyMatchingFunc(const GlobalMuonTrack& mchTrack,
                      const MFTTrack& mftTrack)
{
  auto dx = mchTrack.getX() - mftTrack.getX();
  auto dy = mchTrack.getY() - mftTrack.getY();
  auto score = dx * dx + dy * dy;
  return score;
};

//_________________________________________________________________________________________________
// Sample custom cut criteria that can be passed to MUONMatcher
bool MyMatchingCut(GlobalMuonTrack& mchTrack, MFTTrack& mftTrack)
{
  auto cutDistance = 1.0;
  auto dx = mchTrack.getX() - mftTrack.getX();
  auto dy = mchTrack.getY() - mftTrack.getY();
  auto distance = TMath::Sqrt(dx * dx + dy * dy);
  return distance < cutDistance;
};

//_________________________________________________________________________________________________
// Set Matching function defined by shell variable
void loadAndSetMatchingConfig()
{
  std::string matching_fcn;
  if (gSystem->Getenv("MATCH_SAVE_ALL")) {
    matcher.setMatchSaveAll(true);
    std::cout << " MATCH_SAVE_ALL set to true " << std::endl;
  }

  if (gSystem->Getenv("MATCHING_FCN")) {
    matching_fcn = gSystem->Getenv("MATCHING_FCN");
    std::cout << " MATCHING_FCN: " << matching_fcn << std::endl;

    if (matching_fcn.find("matchXY_") < matching_fcn.length()) {
      std::cout << " Setting " << matching_fcn << std::endl;
      matcher.setMatchingFunction(&MUONMatcher::matchMFT_MCH_TracksXY);
    }
    if (matching_fcn.find("matchXYPhiTanl_") < matching_fcn.length()) {
      std::cout << " Setting " << matching_fcn << std::endl;
      matcher.setMatchingFunction(&MUONMatcher::matchMFT_MCH_TracksXYPhiTanl);
    }
    if (matching_fcn.find("matchALL_") < matching_fcn.length()) {
      std::cout << " Setting " << matching_fcn << std::endl;
      matcher.setMatchingFunction(&MUONMatcher::matchMFT_MCH_TracksAllParam);
    }
    if (matching_fcn.find("Hiroshima_") < matching_fcn.length()) {
      std::cout << " Setting " << matching_fcn << std::endl;
      matcher.setMatchingFunction(&MUONMatcher::Hiroshima);
    }
    if (matching_fcn.find("trainedML_") < matching_fcn.length()) {
      std::cout << " Setting " << matching_fcn << std::endl;
      matcher.setMatchingFunction(&MUONMatcher::matchTrainedML);
      if (gSystem->Getenv("ML_WEIGHTFILE")) {
        float_t scorecut = 0.5;
        if (gSystem->Getenv("ML_SCORECUT"))
          scorecut = atof(gSystem->Getenv("ML_SCORECUT"));
        std::string weightfilename = gSystem->Getenv("ML_WEIGHTFILE");
        matcher.configureTMVA(weightfilename, scorecut);
        std::cout << "Setting TMVA weight file: " << weightfilename
                  << std::endl;
        std::cout << "Setting matching score cut = " << scorecut << std::endl;
      } else {
        std::cout << "Missing TMVA Weight File!" << std::endl;
        exit(1);
      }
    }
  } else {
    std::cout << " Setting default matching function: ";
    matcher.setMatchingFunction(&MUONMatcher::matchMFT_MCH_TracksAllParam);
  }

  if (gSystem->Getenv("MATCHING_PLANEZ")) {
    double matching_planeZ = atof(gSystem->Getenv("MATCHING_PLANEZ"));
    std::cout << " MATCHING_PLANEZ: " << matching_planeZ << std::endl;
    matcher.SetMatchingPlaneZ(matching_planeZ);
  }

  std::string matching_cutfcn;
  if (gSystem->Getenv("MATCHING_CUTFCN")) {
    matching_cutfcn = gSystem->Getenv("MATCHING_CUTFCN");
    std::cout << " MATCHING_CUTFCN: " << matching_cutfcn << std::endl;

    if (matching_cutfcn.find("cutDisabled_") < matching_cutfcn.length()) {
      std::cout << " Setting " << matching_cutfcn << std::endl;
      matcher.setCutFunction(&MUONMatcher::matchCutDisabled);
    }
    if (matching_cutfcn.find("cutDistance_") < matching_cutfcn.length()) {
      std::cout << " Setting " << matching_cutfcn << std::endl;
      matcher.setCutFunction(&MUONMatcher::matchCutDistance);
    }
    if (matching_cutfcn.find("cutDistanceSigma_") < matching_cutfcn.length()) {
      std::cout << " Setting " << matching_cutfcn << std::endl;
      matcher.setCutFunction(&MUONMatcher::matchCutDistanceSigma);
    }
    if (matching_cutfcn.find("cutDistanceAndAngles_") <
        matching_cutfcn.length()) {
      std::cout << " Setting " << matching_cutfcn << std::endl;
      matcher.setCutFunction(&MUONMatcher::matchCutDistanceAndAngles);
    }
    if (matching_cutfcn.find("cutDistanceAndAngles3Sigma_") <
        matching_cutfcn.length()) {
      std::cout << " Setting " << matching_cutfcn << std::endl;
      matcher.setCutFunction(&MUONMatcher::matchCut3SigmaXYAngles);
    }
    if (matching_cutfcn.find("cutDistanceAndAnglesVar_") <
        matching_cutfcn.length()) {
      std::cout << " Setting " << matching_cutfcn << std::endl;
      matcher.setCutFunction(&MUONMatcher::matchCutVarXYAngles);
    }
  }

  if (gSystem->Getenv("ENABLECHARGEMATCHCUT")) {
    matcher.enableChargeMatchCut();
  }

  if (gSystem->Getenv("MATCHING_CUTPARAM0")) {
    double matching_cutparam0 = atof(gSystem->Getenv("MATCHING_CUTPARAM0"));
    std::cout << " MATCHING_CUTPARAM0: " << matching_cutparam0 << std::endl;
    matcher.setCutParam(0, matching_cutparam0);
  }

  if (gSystem->Getenv("MATCHING_CUTPARAM1")) {
    double matching_cutparam1 = atof(gSystem->Getenv("MATCHING_CUTPARAM1"));
    std::cout << " MATCHING_CUTPARAM1: " << matching_cutparam1 << std::endl;
    matcher.setCutParam(1, matching_cutparam1);
  }

  if (gSystem->Getenv("MATCHING_CUTPARAM2")) {
    double matching_cutparam2 = atof(gSystem->Getenv("MATCHING_CUTPARAM2"));
    std::cout << " MATCHING_CUTPARAM2: " << matching_cutparam2 << std::endl;
    matcher.setCutParam(2, matching_cutparam2);
  }

  if (gSystem->Getenv("VERBOSEMATCHING")) {
    std::cout << " Vebose matching enabled." << std::endl;
    matcher.SetVerbosity(true);
  }
}

//_________________________________________________________________________________________________
int evalMLExportOrTrain()
{

  // Runs track matching event-by-event or generate training data
  if (gSystem->Getenv("ML_EXPORTTRAINDATA")) {
    int nMCHTracks = atoi(gSystem->Getenv("ML_EXPORTTRAINDATA"));
    std::cout << " Generate ML traning data file for " << nMCHTracks << " MCH tracks." << std::endl;
    matcher.exportTrainingDataRoot(nMCHTracks);
    exit(0);
  }

  // Runs track matching event-by-event or generate training data
  if (gSystem->Getenv("TRAIN_ML")) {
    matcher.MLTraining();
    exit(0);
  }
  return 0;
}

//_________________________________________________________________________________________________
int runMatching()
{

  // gSystem->Load("libO2MCHTracking");
  // Custom matching function
  // matcher.setCustomMatchingFunction(&MyMatchingFunc,
  // "_aliasForMyMatchingFunction");

  // Define ML input features
  // ML Features defined by lambda expression
  //matcher.setMLFeatureFunction([](const MCHTrackConv& mchTrack, const MFTTrack& mftTrack, float* features) {
  //  features[0] = mftTrack.getX() - mchTrack.getX();
  //  features[1] = mftTrack.getY() - mchTrack.getY();
  //  features[2] = mftTrack.getPhi() - mchTrack.getPhi();
  //  features[3] = mftTrack.getTanl() - mchTrack.getTanl();
  //  features[4] = mftTrack.getInvQPt() - mchTrack.getInvQPt();
  //},
  //                             5, "ML5ParDeltas");

  // ML Features and names defined by separete function(s)
  matcher.setMLFeatureFunction(MLParCov40Features, 40, "MLParCov40Features", MLParCov40FeaturesNames);
  //matcher.setMLFeatureFunction(MLParCovChiNPts42Features, 42, "MLParCovChiNPts42Features", MLParCovChiNPts42FeaturesNames);

  // Configure matcher according command line options
  loadAndSetMatchingConfig();

  // Load MCH tracks
  matcher.loadMCHTracks();

  // Load MFT tracks and propagates to matching plane
  matcher.loadMFTTracksOut();

  // Propagate MCH tracks to matching plane and convert parameters and
  // covariances matrix to MFT coordinate system
  matcher.initGlobalTracks();

  evalMLExportOrTrain();

  matcher.runEventMatching();

  // Kalman filter
  matcher.fitTracks();
  matcher.saveGlobalMuonTracks();

  // Export Matching Plane views in png
  matcher.exportNMatchingPlaneViews(10);

  return 0;
}
