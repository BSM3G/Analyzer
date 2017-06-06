#include "Histo.h"

Histogramer::Histogramer() : outfile(nullptr) {}

Histogramer::Histogramer(int _Npdf, string histname, string cutname, string outfilename, bool _isData, vector<string>& folderCuts): 
  outfile(nullptr), outname(outfilename), Npdf(_Npdf), isData(_isData) {

  read_cuts(cutname, folderCuts);
  NFolders = folders.size();
  read_hist(histname);

  if(folderCuts.size() != 0) {
    CR = true;
    for(auto it: data) it.second->setControlRegions();
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
  CR = rhs.CR;

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
  CR = rhs.CR;
  data.swap(rhs.data);

  rhs.data.clear();
  rhs.outfile = nullptr;

  return *this;
}


Histogramer::Histogramer(const Histogramer& rhs) :
  outname(rhs.outname), NFolders(rhs.NFolders), isData(rhs.isData), CR(rhs.CR)
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
  outname(rhs.outname), NFolders(rhs.NFolders), isData(rhs.isData), CR(rhs.CR)
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
  if(outfile != nullptr)
    outfile->Close();
  
  for(auto it: data_order) {
    delete data[it];
    data[it] = nullptr;
  }
}


void Histogramer::read_hist(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");  

  if(!info_file) {
    cout << "ERROR: Didn't Read Histo File!" << endl;
    cout << filename << endl;
    exit(1);
  }

  vector<string> stemp;
  string group,line;
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
      accept = stoi(stemp[1]) && !(isData && group.find("Gen") != string::npos);
      if(accept) {
	data[group] = new DataBinner();
	data_order.push_back(group); 
      }
    } else if(!accept) continue;
    else if(stemp.size() == 4) {
      string name = extractHistname(group, stemp[0]);
      data[group]->Add_Hist(name, stemp[0], stod(stemp[1]), stod(stemp[2]), stod(stemp[3]), NFolders);
    } else if(stemp.size() == 7) {
      string name = extractHistname(group, stemp[0]);
      data[group]->Add_Hist(name, stemp[0], stod(stemp[1]), stod(stemp[2]), stod(stemp[3]),stod(stemp[4]), stod(stemp[5]), stod(stemp[6]), NFolders);
    }
  }

  info_file.close(); 
}


string Histogramer::extractHistname(string group, string histo) const {
  regex reg ("((Tau|Muon|Electron)+(1|2)+)");
  smatch m;

  string stringkey = group.erase(0,4);
  regex first (stringkey+"(_)?");
  histo = regex_replace(histo,first, "");
  if(stringkey.find("Di") != string::npos) {
    stringkey=stringkey+"1"+stringkey+"2";
  }

  int i=1;
  while(regex_search(stringkey,m,reg)) {
    regex key(m[0].str());
    histo = regex_replace(histo,key,"Part"+to_string(i));
    stringkey = m.suffix().str();
    i++;
  }
  return histo;
}


void Histogramer::read_cuts(string filename, vector<string>& folderCuts) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");  

  if(!info_file) {
    cout << "ERROR: Didn't Read Histo File!" << endl;
    cout << filename << endl;
    exit(1);
  }

  vector<string> stemp;
  string name,line;
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


void Histogramer::fillCRFolderNames(string sofar, int index, bool isFirst, const vector<string>& variables) {
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


void Histogramer::fill_histogram() {

  outfile = new TFile(outname.c_str(), "RECREATE");
  for(auto it: folders) {
    outfile->mkdir( it.c_str() );	
  }
  
  for(auto it: data_order) {
    data[it]->write_histogram(outfile, folders);
  }
  outfile->Close();
}  

void Histogramer::addVal(double value, string group, int maxcut, string histn, double weight) {
  int maxFolder=0;


  if(CR) maxFolder = maxcut;
  else {
    for(int i = 0; i < NFolders; i++) {
      if(maxcut > folderToCutNum[i]) maxFolder++;
      else break;
    }
  }

  data[group]->AddPoint(histn, maxFolder, value, weight);

}

void Histogramer::addVal(double valuex, double valuey, string group, int maxcut, string histn, double weight) {
  int maxFolder=0;

  if(CR) maxFolder = maxcut;
  else {
    for(int i = 0; i < NFolders; i++) {
      if(maxcut > folderToCutNum[i]) maxFolder++;
      else break;
    }
  }

  data[group]->AddPoint(histn, maxFolder, valuex, valuey, weight);
}

    
