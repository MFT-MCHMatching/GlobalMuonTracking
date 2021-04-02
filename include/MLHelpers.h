
//_________________________________________________________________________________________________
string opt_reader(const char *filename = "MLConfigs.xml") {
  // First create engine
  TXMLEngine xml;
  // Now try to parse xml file
  // Only file with restricted xml syntax are supported
  XMLDocPointer_t xmldoc = xml.ParseFile(filename);
  if (!xmldoc)
    exit(1);
  // take access to main node
  XMLNodePointer_t mainnode = xml.DocGetRootElement(xmldoc);
  // display recursively all nodes and subnodes
  // DisplayNode(xml, mainnode, 1);

  std::string training_str("");
  std::string MLLayout;
  std::string MLStrat;
  std::string MLOpt;
  std::string MLMethodType = gSystem->Getenv("TRAIN_ML_METHOD");

  const char *content2;
  bool param_found;
  std::map<string, string> configuration_map;
  std::vector<string> configvect;

  if (gSystem->Getenv("ML_LAYOUT")) {
    MLLayout = gSystem->Getenv("ML_LAYOUT");
    configuration_map[" Network Layout "] = MLLayout;
    //      network_ID += MLLayout + "_";
    configvect.emplace_back(" Network Layout ");
  }
  if (gSystem->Getenv("ML_TRAINING_STRAT")) {
    MLStrat = gSystem->Getenv("ML_TRAINING_STRAT");
    configuration_map[" Training Strategy "] = MLStrat;
    configvect.emplace_back(" Training Strategy ");
    //    network_ID += MLStrat + "_";
  }
  if (gSystem->Getenv("ML_GENERAL_OPT")) {
    MLOpt = gSystem->Getenv("ML_GENERAL_OPT");
    configuration_map[" General Option "] = MLOpt;
    configvect.emplace_back(" General Option ");
    //    network_ID += MLOpt + "_";
  }

  while (xml.GetNodeName(mainnode) != MLMethodType) {
    mainnode = xml.GetNext(mainnode);
  }
  for (auto config_str : configvect) {

    param_found = false;
    XMLNodePointer_t config =
        xml.GetChild(mainnode); // Note that the program will search the entire
                                // .xml file for the entered option
    while (config) {
      XMLNodePointer_t param =
          xml.GetChild(config); // second child layers should be the parameters
      while (param) {
        if (configuration_map[config_str] == xml.GetNodeName(param)) {
          content2 = xml.GetNodeContent(param);
          TString paramstr(content2);
          std::cout << config_str << xml.GetNodeName(param) << " found! "
                    << std::endl;
          training_str += paramstr + ":";
          param_found = true;
          break;
        }
        param = xml.GetNext(param);
      }
      if (param_found)
        break;
      config = xml.GetNext(config);
    }
    if (!param_found) {
      std::cout << config_str << configuration_map[config_str]
                << " NOT found! Using TMVA's default where needed. \n"
                << std::endl;
    }
  }
  // Release memory before exit
  xml.FreeDoc(xmldoc);

  return training_str;
}
