#include "Histo.h"
#include "unistd.h"
#include "Compression.h"

Histogramer::Histogramer() : outfile(nullptr) {}

Histogramer::Histogramer(int _Npdf, std::string histname, std::string cutname, std::string outfilename, bool _isData, std::vector<std::string>& folderCuts, const std::vector<std::string> &syst_unvertainties ): outname(outfilename), 
outfile(nullptr), Npdf(_Npdf), isData(_isData) {

  //no syst uncertainty hist object
  if (syst_unvertainties.size()==0){
    read_cuts(cutname, folderCuts);
    //outfile = new TFile(outfilename.c_str(), "RECREATE");
  }else{
    read_syst(syst_unvertainties);
  }
  NFolders = folders.size();
  read_hist(histname);

  if(folderCuts.size() != 0 || syst_unvertainties.size() != 0) {
    fillSingle = true;
    for(auto it: data) it.second->setSingleFill();
  }
}


Histogramer& Histogramer::operator=(const Histogramer& rhs) {
  if(this == &rhs) return *this;

  outname = rhs.outname;
  NFolders = rhs.NFolders;
  isData = rhs.isData;

  cuts = rhs.cuts;
  cut_order = rhs.cut_order;
  folders = rhs.folders;
  folderToCutNum = rhs.folderToCutNum;
  data_order.reserve(rhs.data_order.size());
  data_order = rhs.data_order;
  fillSingle = rhs.fillSingle;

  for(auto mit: rhs.data) {
    data[mit.first] = new DataBinner(*(mit.second));
  }
  if(rhs.outfile != nullptr) {
    outfile = (TFile*)rhs.outfile->Clone();
  }

  return *this;
}

Histogramer& Histogramer::operator=(Histogramer&& rhs) {
  if(this == &rhs) return *this;

  outname = rhs.outname;
  NFolders = rhs.NFolders;
  isData = rhs.isData;
  outfile = rhs.outfile;

  cuts = rhs.cuts;
  cut_order = rhs.cut_order;
  folders = rhs.folders;
  folderToCutNum = rhs.folderToCutNum;

  data_order = rhs.data_order;
  fillSingle = rhs.fillSingle;
  data.swap(rhs.data);

  rhs.data.clear();
  rhs.outfile = nullptr;

  return *this;
}


Histogramer::Histogramer(const Histogramer& rhs) :
outname(rhs.outname), NFolders(rhs.NFolders), isData(rhs.isData), fillSingle(rhs.fillSingle)
{
  cuts = rhs.cuts;
  cut_order = rhs.cut_order;
  folders = rhs.folders;
  folderToCutNum = rhs.folderToCutNum;
  data_order = rhs.data_order;

  for(auto mit: rhs.data) {
    data[mit.first] = new DataBinner(*(mit.second));
  }
  if(rhs.outfile != nullptr) {
    outfile = (TFile*)rhs.outfile->Clone();
  }

}

Histogramer::Histogramer(Histogramer&& rhs) :
outname(rhs.outname), NFolders(rhs.NFolders), isData(rhs.isData), fillSingle(rhs.fillSingle)
{
  cuts = rhs.cuts;
  cut_order = rhs.cut_order;
  folders = rhs.folders;
  folderToCutNum = rhs.folderToCutNum;
  data_order = rhs.data_order;
  data.swap(rhs.data);
  outfile = rhs.outfile;

  rhs.data.clear();
  outfile = nullptr;
}


Histogramer::~Histogramer() {
  //if(outfile != nullptr)
    //outfile->Close();
}


void Histogramer::read_hist(std::string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  std::ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "ERROR: Didn't Read Histo File!" << std::endl;
    std::cout << filename << std::endl;
    exit(1);
  }

  std::vector<std::string> stemp;
  std::string group,line;
  bool accept = false;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }

    if(stemp.size() == 0) continue;
    else if(stemp.size() == 2) {
      group = stemp[0];
      accept = stoi(stemp[1]) && !(isData && group.find("Gen") != std::string::npos);
      if(accept) {
        data[group] = new DataBinner();
        data_order.push_back(group);
      }
    } else if(!accept) continue;
    else if(stemp.size() == 4) {
      std::string name = extractHistname(group, stemp[0]);
      data[group]->Add_Hist(name, stemp[0], stod(stemp[1]), stod(stemp[2]), stod(stemp[3]), NFolders);
    } else if(stemp.size() == 7) {
      std::string name = extractHistname(group, stemp[0]);
      data[group]->Add_Hist(name, stemp[0], stod(stemp[1]), stod(stemp[2]), stod(stemp[3]),stod(stemp[4]), stod(stemp[5]), stod(stemp[6]), NFolders);
    }
  }

  info_file.close();
  
  data["Eff"] = new DataBinner();
  std::vector<std::string> allLepNames={"Electron","Muon","Tau"};
  for(std::string s : allLepNames){
    std::string name="eff_"+s;
    data["Eff"]->Add_Hist(name+"Pt",  300, 0, 3000, 1);
    data["Eff"]->Add_Hist(name+"Eta", 100, -5, 5, 1);
    data["Eff"]->Add_Hist(name+"Phi", 100, -3.14159, 3.14159, 1);
    name="eff_Reco_"+s;
    data["Eff"]->Add_Hist(name+"Pt",  300, 0, 3000, 1);
    data["Eff"]->Add_Hist(name+"Eta", 100, -5, 5, 1);
    data["Eff"]->Add_Hist(name+"Phi", 100, -3.14159, 3.14159, 1);
  }
  data_order.push_back("Eff");
}


std::string Histogramer::extractHistname(std::string group, std::string histo) const {
  std::regex reg ("((Tau|Muon|Electron)+(1|2)+)");
  std::smatch m;

  std::string stringkey = group.erase(0,4);
  std::regex first (stringkey+"(_)?");
  histo = std::regex_replace(histo,first, "");
  if(stringkey.find("Di") != std::string::npos) {
    stringkey=stringkey+"1"+stringkey+"2";
  }

  int i=1;
  while(std::regex_search(stringkey,m,reg)) {
    std::regex key(m[0].str());
    histo = std::regex_replace(histo,key,"Part"+std::to_string(i));
    stringkey = m.suffix().str();
    i++;
  }
  return histo;
}


void Histogramer::read_cuts(std::string filename, std::vector<std::string>& folderCuts) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  std::ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "ERROR: Didn't Read Histo File!" << std::endl;
    std::cout << filename << std::endl;
    exit(1);
  }

  std::vector<std::string> stemp;
  std::string name,line;
  int i = 0;

  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }

    if(stemp.size() == 3) {
      name = stemp[0];
      if(name[0]=='*' && name[1]=='*' && name[2]=='*') {
        name.erase(0,3);
        if(folderCuts.size() == 0) {
          folders.push_back(name);
          folderToCutNum.push_back(i);
        }
      } else if(stemp[1] == "0" && stemp[2] == "-1") continue;   ////remove unnecessary cuts
      cuts[name] = std::make_pair(stoi(stemp[1]),stoi(stemp[2]));
      cut_order.push_back(name);
      i++;
    }
  }

  if(folderCuts.size() != 0) {
    fillCRFolderNames("", 0, true, folderCuts);
  } else if(folders.size() == 0 && cut_order.size() == 0){
    folders.push_back(name);
    folderToCutNum.push_back(i-1);
  } else if(cut_order.size() != 0 && ( folders.size() == 0 || folders.back() != cut_order.back())) {
    folders.push_back(cut_order.back());
    folderToCutNum.push_back(cut_order.size() -1);
  }

  info_file.close();
}


void Histogramer::read_syst(const std::vector<std::string>& syst_uncertainties) {

  int i=0;
  for(const std::string &syst : syst_uncertainties){
    folders.push_back(syst);
    folderToCutNum.push_back(i);
    i++;
  }

}


void Histogramer::fillCRFolderNames(std::string sofar, int index, bool isFirst, const std::vector<std::string>& variables) {
  if(index >= (int)variables.size()) {
    folders.push_back(sofar);
    return;
  }
  if(isFirst) {
    fillCRFolderNames(sofar+variables[index]+"<"+variables[index+1]+"_", index+2, true, variables);
    fillCRFolderNames(sofar, index, false, variables);
  } else {
    fillCRFolderNames(sofar+variables[index]+">"+variables[index+1]+"_", index+2, true, variables);
  }
}


void Histogramer::fill_histogram(TFile* _outfile, std::string subfolder) {
  //if( access( outname.c_str(), F_OK ) == -1 ){
    //outfile = new TFile(outname.c_str(), "RECREATE", outname.c_str(), ROOT::CompressionSettings(ROOT::kLZMA, 9));
  //}else{
    //outfile = new TFile(outname.c_str(), "UPDATE", outname.c_str(), ROOT::CompressionSettings(ROOT::kLZMA, 9));
  //}
  outfile = _outfile;
  if(subfolder!=""){
    outfile->mkdir(subfolder.c_str());
    outfile->cd(subfolder.c_str());
    for(std::string it: folders) {
      outfile->mkdir( (subfolder+"/"+it).c_str() );
    }
  }else{
    for(auto it: folders) {
      outfile->mkdir( it.c_str() );
    }
    if(outfile->GetDirectory("Eff")==nullptr)
    outfile->mkdir("Eff");
  }
  for(auto it: data_order) {
    data[it]->write_histogram(outfile, folders, subfolder);
  }
  outfile->cd();
  if(subfolder!=""){
    outfile->cd(subfolder.c_str());
  }
  for (std::unordered_map<std::string, TTree * >::iterator it = trees.begin(); it != trees.end(); ++it) {
    it->second->Write();
  }
  //outfile->Close();
}

void Histogramer::createTree(std::unordered_map< std::string , float > *m, std::string name){
  trees[name] = new TTree(name.c_str(), name.c_str());
  for (std::unordered_map< std::string , float >::iterator it = m->begin(); it != m->end(); it++) {
    trees[name]->Branch(it->first.c_str(), &(it->second), (it->first+"/F").c_str());
  }
}

void Histogramer::fillTree(std::string name) {
  trees[name]->Fill();
}


void Histogramer::addVal(double valuex, double valuey, std::string group, int maxcut, std::string histn, double weight) {
  int maxFolder=0;


  if(fillSingle) maxFolder = maxcut;
  else {
    for(int i = 0; i < NFolders; i++) {
      if(maxcut > folderToCutNum[i]) maxFolder++;
      else break;
    }
  }
  data[group]->AddPoint(histn, maxFolder, valuex, valuey, weight);
}

void Histogramer::addVal(double value, std::string group, int maxcut, std::string histn, double weight) {
  int maxFolder=0;


  if(fillSingle) maxFolder = maxcut;
  else {
    for(int i = 0; i < NFolders; i++) {
      if(maxcut > folderToCutNum[i]) maxFolder++;
      else break;
    }
  }
  data[group]->AddPoint(histn, maxFolder, value, weight);
}



void Histogramer::addEffiency(std::string histn ,double value ,bool passFail,int maxFolder=0){
  
  data["Eff"]->AddEff(histn, maxFolder, value,passFail);
}


