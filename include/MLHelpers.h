
//_________________________________________________________________________________________________
string DNN_read(const std::string& MLLayout, const std::string& MLStrat, const std::string& MLOpt, const char* filename = "MLConfigs.xml")
{
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

  std::map<string, string> configuration_map;
  configuration_map[" Network Layout "] = MLLayout;
  configuration_map[" Training Strategy "] = MLStrat;
  configuration_map[" General Option "] = MLOpt;

  const char* content2;
  XMLNodePointer_t config = xml.GetChild(mainnode); //first config (first child) should be layouts, second training strat and so on

  std::vector<string> configvect{" Network Layout ", " Training Strategy ", " General Option "};
  for (auto config_str : configvect) {
    XMLNodePointer_t param = xml.GetChild(config); //second child layers should be the parameters
    while (param) {
      if (configuration_map[config_str] == xml.GetNodeName(param)) {
        content2 = xml.GetNodeContent(param);
        TString paramstr(content2);
        std::cout << config_str << xml.GetNodeName(param) << " found! " << std::endl;
        //        std::cout<<" Content: " << paramstr<<"\n"<<std::endl;
        training_str += paramstr + ":";
        break;
      }
      param = xml.GetNext(param);
      if (!param) {
        std::cout << config_str << configuration_map[config_str] << " NOT found! Aborting...\n"
                  << std::endl;
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
