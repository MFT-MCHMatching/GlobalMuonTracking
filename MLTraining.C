
/// Macro for training a regression DNN (DL) for MCH-MFT Track Matching
/// It was made based on the example macro of TMVA's Regression, authored by Andreas Hoecker
///								https://root.cern/doc/master/TMVARegression_8C.html



#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"

#include "TXMLEngine.h"

using namespace TMVA;

string DNN_read(const std::string &MLLayout, const std::string &MLStrat, const std::string &MLOpt, const char* filename = "MLConfigs.xml")
{
   // First create engine
   TXMLEngine xml;
   // Now try to parse xml file
   // Only file with restricted xml syntax are supported
   XMLDocPointer_t xmldoc = xml.ParseFile(filename);
   if (!xmldoc) exit(1);
   // take access to main node
   XMLNodePointer_t mainnode = xml.DocGetRootElement(xmldoc);
   // display recursively all nodes and subnodes
  // DisplayNode(xml, mainnode, 1);

  std::string training_str("");

  std::map<string, string> configuration_map;
  configuration_map[" Network Layout "] = MLLayout;
  configuration_map[" Training Strategy "] = MLStrat;
  configuration_map[" General Option "] = MLOpt;

  const char *content2;
  XMLNodePointer_t config = xml.GetChild(mainnode);  //first config (first child) should be layouts, second training strat and so on

  std::vector<string> configvect{" Network Layout "," Training Strategy "," General Option "};
  for (auto config_str:configvect) {
    XMLNodePointer_t param = xml.GetChild(config);			//second child layers should be the parameters
    while(param) {
      if(configuration_map[config_str] == xml.GetNodeName(param)) {
        content2 = xml.GetNodeContent(param);
        TString paramstr(content2);
        std::cout<< config_str << xml.GetNodeName(param)<<" found! "<<std::endl;
//        std::cout<<" Content: " << paramstr<<"\n"<<std::endl;
        training_str += paramstr+":";
        break;
      }
      param = xml.GetNext(param);
      if(!param) {
        std::cout << config_str << configuration_map[config_str] << " NOT found! Aborting...\n" << std::endl;
        exit(1);
      }
	}	
    config = xml.GetNext(config);
  }
//		std::cout<<" Complete ML String: "<<training_str<<"\n"<<std::endl;
   // Release memory before exit
  xml.FreeDoc(xmldoc);

return training_str;
}

void DLRegression(std::string input_name, std::string trainingfile, std::string trainingstr)
{
  
	//---------------------------------------------------------------
	// This loads the library
  TMVA::Tools::Instance();

  std::cout << "==> Start TMVARegression" << std::endl;

	TFile *input(0);
	if (!gSystem->AccessPathName( trainingfile.c_str() )) {    //TODO pass file check to shell script. Must solve problem with name cuts: ln128 -131
		input = TFile::Open( trainingfile.c_str() );
	}
	if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
  }
  std::cout << "--- TMVARegression           : Using training input file: " << input->GetName() << std::endl;

  // Register the regression tree
  TTree *regTree = (TTree*)input->Get("matchTree");

  // Let's initialize the factory object (analysis class)...:
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results.
  // The weightfile name will be the the combined strings of the two arguments
  //
  // In the third argument we can set the analysis type: classification or regression

	TString methodname = input_name.c_str();

  std::size_t pos1 = trainingfile.find("MLTraining_");
  std::size_t pos = trainingfile.find(".root");
  methodname += "_" + trainingfile.substr(pos1,pos);
  TString outfileName( methodname+".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  TMVA::Factory* factory = new TMVA::Factory( "Trained_ML", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
				

  // Dataloader object - this will handle the data (The argument also defines the name of the directory containing the weights' file)
   TMVA::DataLoader* dataloader=new TMVA::DataLoader("trainedMLs");

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
	// First argument Must be the same as in the TTree. In the last argument the appropriate type must be set 
  dataloader->AddVariable( "MFT_X", "MFT_X_pos", "units", 'F' );
  dataloader->AddVariable( "MFT_Y", "MFT_Y_pos", "units", 'F' );
  dataloader->AddVariable( "MFT_Phi", "MFT_Phi", "units", 'F' );
  dataloader->AddVariable( "MFT_Tanl", "MFT_Tanl", "units", 'F' );
  dataloader->AddVariable( "MFT_InvQPt", "MFT_InvQPt", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov00", "MFT_Cov00", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov01", "MFT_Cov01", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov11", "MFT_Cov11", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov02", "MFT_Cov02", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov12", "MFT_Cov12", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov22", "MFT_Cov22", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov03", "MFT_Cov03", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov13", "MFT_Cov13", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov23", "MFT_Cov23", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov33", "MFT_Cov33", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov04", "MFT_Cov04", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov14", "MFT_Cov14", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov24", "MFT_Cov24", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov34", "MFT_Cov34", "units", 'F' );
  dataloader->AddVariable( "MFT_Cov44", "MFT_Cov44", "units", 'F' );
  dataloader->AddVariable( "MCH_X", "MCH_X_pos", "units", 'F' );
  dataloader->AddVariable( "MCH_Y", "MCH_Y_pos", "units", 'F' );
  dataloader->AddVariable( "MCH_Phi", "MCH_Phi", "units", 'F' );
  dataloader->AddVariable( "MCH_Tanl", "MCH_Tanl", "units", 'F' );
  dataloader->AddVariable( "MCH_InvQPt", "MCH_InvQPt", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov00", "MCH_Cov00", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov01", "MCH_Cov01", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov11", "MCH_Cov11", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov02", "MCH_Cov02", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov12", "MCH_Cov12", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov22", "MCH_Cov22", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov03", "MCH_Cov03", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov13", "MCH_Cov13", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov23", "MCH_Cov23", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov33", "MCH_Cov33", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov04", "MCH_Cov04", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov14", "MCH_Cov14", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov24", "MCH_Cov24", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov34", "MCH_Cov34", "units", 'F' );
  dataloader->AddVariable( "MCH_Cov44", "MCH_Cov44", "units", 'F' );
   // Add the variable carrying the regression target
  dataloader->AddTarget( "Truth" );
 
   // global event weights per tree (see below for setting event-wise weights)
  Double_t regWeight  = 1.0;

   // You can now add the tree to the dataloader:
  dataloader->AddRegressionTree( regTree, regWeight, Types::kTraining );
   // ** In case we need some data for test, we shoudl use:
//  dataloader->AddRegressionTree( testTree, regWeight, Types::kTesting);

	// Last, but not least, lets set the options for training (and testing):

   // Apply additional cuts on the signal and background samples (can be different)
  TCut mycut = ""; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";

  dataloader->PrepareTrainingAndTestTree( mycut,
                                         "SplitMode=Random:NormMode=NumEvents:!V" );

	// Book (SAVE) the method
  factory->BookMethod(dataloader, TMVA::Types::kDL, methodname, trainingstr); 


   // --------------------------------------------------------------------------------------------------
  std::ofstream info("time_"+methodname+".txt");
  TStopwatch *timewatch = new TStopwatch();
   // Train MVAs using the set of training events
  factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
//   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
//   factory->EvaluateAllMethods();

  info<<"\n END OF THE TRAINING \n";
  info<< " ~CPU time (s): " << timewatch->CpuTime()<<"\n\n"<< " ~Real time(s): "<< timewatch->RealTime()<<"\n\n";
   // --------------------------------------------------------------

   // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVARegression is done!" << std::endl;

  delete factory;
  delete dataloader;

   // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );
}

void MLTraining()
{
  std::string MLLayout = gSystem->Getenv("ML_LAYOUT");
  std::string MLStrat = gSystem->Getenv("ML_TRAINING_STRAT");
  std::string MLOpt = gSystem->Getenv("ML_GENERAL_OPT");
  std::string training_file = gSystem->Getenv("ML_TRAINING_FILE");
	
  std::string training_string ( DNN_read(MLLayout, MLStrat, MLOpt) );
  std::string network_ID (MLLayout + "_" + MLStrat + "_" + MLOpt);
  std::cout<<" Network name: "<< network_ID<< "\n"<<std::endl;
//	std::cout<<" Training file "<< training_file<< "\n"<<std::endl;
  DLRegression(network_ID,training_file,training_string);
}



