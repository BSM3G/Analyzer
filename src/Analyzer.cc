#include "Analyzer.h"
//// Used to convert Enums to integers
#define ival(x) static_cast<int>(x)
//// BIG_NUM = sqrt(sizeof(int)) so can use diparticle convention of
//// index = BIG_NUM * i1 + i2
//// This insures easy way to extract indices
//// Needs to be changed if go to size_t instead (if want to play safe
#define BIG_NUM 46340

///// Macros defined to shorten code.  Made since lines used A LOT and repeative.  May change to inlines
///// if tests show no loss in speed
#define histAddVal2(val1, val2, name) histo.addVal(val1, val2, group, max, name, wgt)
#define histAddVal(val, name) histo.addVal(val, group, max, name, wgt)
#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

typedef vector<int>::iterator vec_iter;

//////////////////////////////////////////////////////////////////
///////////////////CONSTANTS DEFINITONS///////////////////////////
//////////////////////////////////////////////////////////////////

//Filespace that has all of the .in files
const string FILESPACE = "PartDet/";
const string PUSPACE = "Pileup/";

const vector<CUTS> Analyzer::genCuts = {
  CUTS::eGTau, CUTS::eNuTau, CUTS::eGTop, 
  CUTS::eGElec, CUTS::eGMuon, CUTS::eGZ, 
  CUTS::eGW, CUTS::eGHiggs
};   

const vector<CUTS> Analyzer::jetCuts = {
  CUTS::eRJet1,  CUTS::eRJet2,   CUTS::eRCenJet,
  CUTS::eR1stJet, CUTS::eR2ndJet, CUTS::eRBJet
};   




const unordered_map<CUTS, vector<CUTS>, EnumHash> Analyzer::adjList = {
  {CUTS::eMuon1Tau1, {CUTS::eRMuon1, CUTS::eRTau1}},
  {CUTS::eMuon1Tau2, {CUTS::eRMuon1, CUTS::eRTau2}},
  {CUTS::eMuon2Tau1, {CUTS::eRMuon2, CUTS::eRTau1}},
  {CUTS::eMuon2Tau2, {CUTS::eRMuon2, CUTS::eRTau2}},

  {CUTS::eElec1Tau1, {CUTS::eRElec1, CUTS::eRTau1}},
  {CUTS::eElec1Tau2, {CUTS::eRElec1, CUTS::eRTau2}},
  {CUTS::eElec2Tau1, {CUTS::eRElec2, CUTS::eRTau1}},
  {CUTS::eMuon2Tau2, {CUTS::eRElec2, CUTS::eRTau2}},

  {CUTS::eMuon1Elec1, {CUTS::eRMuon1, CUTS::eRElec1}},
  {CUTS::eMuon1Elec2, {CUTS::eRMuon1, CUTS::eRElec2}},
  {CUTS::eMuon2Elec1, {CUTS::eRMuon2, CUTS::eRElec1}},
  {CUTS::eMuon2Elec2, {CUTS::eRMuon2, CUTS::eRElec2}},

  {CUTS::eDiElec, {CUTS::eRElec1, CUTS::eRElec2}},
  {CUTS::eDiMuon, {CUTS::eRMuon1, CUTS::eRMuon2, CUTS::eR1stJet, CUTS::eR2ndJet}},
  {CUTS::eDiTau, {CUTS::eRTau1, CUTS::eRTau2}},
  {CUTS::eDiJet, {CUTS::eRJet1, CUTS::eRJet2}},
  {CUTS::eSusyCom, {CUTS::eR1stJet, CUTS::eR2ndJet}},
  {CUTS::eGTau, {CUTS::eNuTau}},
  {CUTS::eGen, genCuts}
};


const unordered_map<string, CUTS> Analyzer::cut_num = { 
  {"NGenTau", CUTS::eGTau},                             {"NGenTop", CUTS::eGTop}, 
  {"NGenElectron", CUTS::eGElec},                       {"NGenMuon", CUTS::eGMuon}, 
  {"NGenZ", CUTS::eGZ},                                 {"NGenW", CUTS::eGW}, 
  {"NGenHiggs", CUTS::eGHiggs},                         {"NRecoVertex", CUTS::eRVertex}, 
  {"NRecoMuon1", CUTS::eRMuon1},                        {"NRecoMuon2", CUTS::eRMuon2}, 
  {"NRecoElectron1", CUTS::eRElec1},                    {"NRecoElectron2",CUTS::eRElec2},
  {"NRecoTau1", CUTS::eRTau1},                          {"NRecoTau2", CUTS::eRTau2}, 
  {"NRecoJet1", CUTS::eRJet1},                          {"NRecoJet2", CUTS::eRJet2}, 
  {"NRecoCentralJet", CUTS::eRCenJet},                  {"NRecoBJet", CUTS::eRBJet}, 
  {"NRecoTriggers1", CUTS::eRTrig1},                    {"NRecoTriggers2", CUTS::eRTrig2}, 
  {"NRecoFirstLeadingJet", CUTS::eR1stJet},             {"NRecoSecondLeadingJet", CUTS::eR2ndJet},
  {"NDiMuonCombinations", CUTS::eDiMuon},               {"NDiElectronCombinations", CUTS::eDiElec}, 
  {"NDiTauCombinations", CUTS::eDiTau},                 {"NDiJetCombinations", CUTS::eDiJet},
  {"NMuon1Tau1Combinations", CUTS::eMuon1Tau1},         {"NMuon1Tau2Combinations", CUTS::eMuon1Tau2}, 
  {"NMuon2Tau1Combinations", CUTS::eMuon2Tau1},         {"NMuon2Tau2Combinations", CUTS::eMuon2Tau2},
  {"NElectron1Tau1Combinations", CUTS::eElec1Tau1},     {"NElectron1Tau2Combinations", CUTS::eElec1Tau2},
  {"NElectron2Tau1Combinations", CUTS::eElec2Tau1},     {"NElectron2Tau2Combinations", CUTS::eElec2Tau2},
  {"NMuon1Electron1Combinations", CUTS::eMuon1Elec1},   {"NMuon1Electron2Combinations", CUTS::eMuon1Elec2},
  {"NMuon2Electron1Combinations", CUTS::eMuon2Elec1},   {"NMuon2Electron2Combinations", CUTS::eMuon2Elec2},
  {"NLeadJetCombinations", CUTS::eSusyCom},             {"METCut", CUTS::eMET} 
};


//////////////////////////////////////////////////////
//////////////////PUBLIC FUNCTIONS////////////////////
//////////////////////////////////////////////////////

///Constructor
Analyzer::Analyzer(string infile, string outfile, bool setCR) : goodParts(getArray()) {
  cout << "setup start" << endl;
  f = TFile::Open(infile.c_str());
  f->cd("TNT");
  BOOM = (TTree*)f->Get("TNT/BOOM");
  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "TOTAL EVENTS: " << nentries << std::endl;

  for(int i=0; i < nTrigReq; i++) {
    vector<int>* tmpi = new vector<int>();
    vector<string>* tmps = new vector<string>();
    trigPlace[i] = tmpi;
    trigName[i] = tmps;
  }

  setupGeneral(BOOM,infile);

  isData = distats["Run"].bmap.at("isData");
  CalculatePUSystematics = distats["Run"].bmap.at("CalculatePUSystematics");
  initializePileupInfo(distats["Run"].smap.at("MCHistos"), distats["Run"].smap.at("DataHistos"));


  if(!isData) {
    _Gen = new Generated(BOOM, FILESPACE + "Gen_info.in");
  }
  _Electron = new Electron(BOOM, FILESPACE + "Electron_info.in");
  _Muon = new Muon(BOOM, FILESPACE + "Muon_info.in");
  _Tau = new Taus(BOOM, FILESPACE + "Tau_info.in");
  _Jet = new Jet(BOOM, FILESPACE + "Jet_info.in");

  _Electron->findExtraCuts();
  _Muon->findExtraCuts();
  _Tau->findExtraCuts();
  _Jet->findExtraCuts();


  vector<string> cr_variables;
  if(setCR) {
    char buf[64];
    read_info(FILESPACE + "Control_Regions.in");
    crbins = pow(2.0, distats["Control_Region"].dmap.size());
    for(auto maper: distats["Control_Region"].dmap) {
      cr_variables.push_back(maper.first);
      sprintf(buf, "%.*G", 16, maper.second);
      cr_variables.push_back(buf);
    }
    if(isData) {
      if(distats["Control_Region"].smap.find("SR") == distats["Control_Region"].smap.end()) {
	cout << "Using Control Regions with data, but no signal region specified can lead to accidentially unblinding a study  before it should be.  Please specify a SR in the file PartDet/Control_Region.in" << endl;
	exit(1);
      } else if(distats["Control_Region"].smap.at("SR").length() != distats["Control_Region"].dmap.size()) {
	cout << "Signal Region specified incorrectly: check signal region variable to make sure the number of variables matches the number of signs in SR" << endl;
	exit(1);
      }
      int factor = 1;
      SignalRegion = 0;
      for(auto gtltSign: distats["Control_Region"].smap["SR"]) {
	if(gtltSign == '>') SignalRegion += factor;
	factor *= 2;
      }
      if(distats["Control_Region"].bmap.find("Unblind") != distats["Control_Region"].bmap.end()) {
	blinded = !distats["Control_Region"].bmap["Unblind"];
      }
    }
  }

  histo = Histogramer(1, FILESPACE+"Hist_entries.in", FILESPACE+"Cuts.in", outfile, isData, cr_variables);

  if(setCR) {
    cuts_per.resize(histo.get_folders()->size());
    cuts_cumul.resize(histo.get_folders()->size());
  } else {
    cuts_per.resize(histo.get_cuts()->size());
    cuts_cumul.resize(histo.get_cuts()->size());
  }
  create_fillInfo();
  for(auto maper: distats["Control_Region"].dmap) {

    setupCR(maper.first, maper.second);
  }

  setCutNeeds();  
  //  exit(1);
  std::cout << "setup complete" << std::endl << endl;
}

unordered_map<CUTS, vector<int>*, EnumHash> Analyzer::getArray() {
  unordered_map<CUTS, vector<int>*, EnumHash> rmap;
  for(auto e: Enum<CUTS>()) {
    rmap[e] = new vector<int>();
  }
  return rmap;
}



void Analyzer::create_fillInfo() {

  fillInfo["FillLeadingJet"] =    new FillVals(CUTS::eSusyCom, FILLER::Dipart, _Jet, _Jet);
  fillInfo["FillGen"] =   new FillVals(CUTS::eGen, FILLER::Single, _Gen);
  fillInfo["FillTau1"] =    new FillVals(CUTS::eRTau1, FILLER::Single, _Tau);
  fillInfo["FillTau2"] =    new FillVals(CUTS::eRTau2, FILLER::Single, _Tau);
  fillInfo["FillMuon1"] =      new FillVals(CUTS::eRMuon1, FILLER::Single, _Muon);
  fillInfo["FillMuon2"] =      new FillVals(CUTS::eRMuon2, FILLER::Single, _Muon);
  fillInfo["FillElectron1"] =  new FillVals(CUTS::eRElec1, FILLER::Single, _Electron);
  fillInfo["FillElectron2"] =  new FillVals(CUTS::eRElec2, FILLER::Single, _Electron);

  fillInfo["FillJet1"] =       new FillVals(CUTS::eRJet1, FILLER::Single, _Jet);
  fillInfo["FillJet2"] =       new FillVals(CUTS::eRJet2, FILLER::Single, _Jet);
  fillInfo["FillBJet"] =       new FillVals(CUTS::eRBJet, FILLER::Single, _Jet);
  fillInfo["FillCentralJet"] = new FillVals(CUTS::eRCenJet, FILLER::Single, _Jet);

  fillInfo["FillDiMuon"] =     new FillVals(CUTS::eDiMuon, FILLER::Dipart, _Muon, _Muon);
  fillInfo["FillDiTau"] =      new FillVals(CUTS::eDiTau, FILLER::Dipart, _Tau, _Tau);
  fillInfo["FillMetCuts"] =  new FillVals();
  fillInfo["FillDiJet"] =  new FillVals(CUTS::eDiJet, FILLER::Dipart, _Jet, _Jet);

  fillInfo["FillMuon1Tau1"] =  new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon1Tau2"] =  new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau1"] =  new FillVals(CUTS::eMuon2Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau2"] =  new FillVals(CUTS::eMuon2Tau2, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillElectron1Tau1"] =  new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron1Tau2"] =  new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau1"] =  new FillVals(CUTS::eElec2Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau2"] =  new FillVals(CUTS::eElec2Tau2, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillMuon1Electron1"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon1Electron2"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron1"] =  new FillVals(CUTS::eMuon2Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron2"] =  new FillVals(CUTS::eMuon2Elec2, FILLER::Dipart, _Muon, _Electron);

  //////I hate this solution so much.  Its terrible
  fillInfo["FillElectron1Electron2"] =     new FillVals(CUTS::eDiElec, FILLER::Single, _Electron, _Electron);
  fillInfo["FillMuon1Muon2"] =     new FillVals(CUTS::eDiMuon, FILLER::Single, _Muon, _Muon);
  fillInfo["FillTau1Tau2"] =      new FillVals(CUTS::eDiTau, FILLER::Single, _Tau, _Tau);
  
  

  for(auto it: *histo.get_groups()) {
    if(fillInfo[it] == nullptr) fillInfo[it] = new FillVals();
  }

  

}

void Analyzer::setupCR(string var, double val) {
  smatch m;
  regex part ("^(.+)_(.+)$");
  if(regex_match(var, m, part)) {
    string name = m[1];
    string cut = "Fill" + name;
    if(fillInfo.find(cut) == fillInfo.end()) {
      cout << cut << " not found, put into fillInfo" << endl;
      exit(1);
    }
    cout << cut << " " << m[2] << " " << val << " " << name << endl;
    testVec.push_back(new CRTester(fillInfo.at(cut), m[2], val, name));
  } else {
    cout << "Could not process line: " << var << endl;
    exit(1);
  } 

}




////destructor
Analyzer::~Analyzer() {
  delete f;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
  if(!isData) delete _Gen;

  for(auto pair: fillInfo) {
    delete pair.second;
  }
  
  for(auto e: Enum<CUTS>()) {
    delete goodParts[e];
  }

  for(int i=0; i < nTrigReq; i++) {
    delete trigPlace[i];
    delete trigName[i];
  }
}


///resets values so analysis can start
void Analyzer::clear_values() {
  for(auto e: Enum<CUTS>()) {
    goodParts[e]->clear();
  }
  deltaMEx=0;
  deltaMEy=0;
  sumpxForMht = 0.0;
  sumpyForMht = 0.0;
  sumptForHt  = 0.0;
  leadIndex=-1;
  maxCut = 0;
}

///Function that does most of the work.  Calculates the number of each particle
void Analyzer::preprocess(int event) {
  BOOM->GetEntry(event);
  theMETVector.SetPxPyPzE(Met[0], Met[1], Met[2], sqrt(pow(Met[0],2) + pow(Met[1],2)));
  pu_weight = (!isData && CalculatePUSystematics) ? hPU[(int)(nTruePU+1)] : 1.0;

  // SET NUMBER OF GEN PARTICLES
  if(!isData){
    getGoodGen(_Gen->pstats["Gen"]);
    getGoodTauNu();
  }

  //////Smearing  
  smearLepton(*_Electron, CUTS::eGElec, _Electron->pstats["Smear"]);
  smearLepton(*_Muon, CUTS::eGMuon, _Muon->pstats["Smear"]);
  smearLepton(*_Tau, CUTS::eGTau, _Tau->pstats["Smear"]);
  smearJet(_Jet->pstats["Smear"]);
  
  //////Triggers and Vertices
  goodParts[CUTS::eRVertex]->resize(bestVertices);
  TriggerCuts(*(trigPlace[0]), *(trigName[0]), CUTS::eRTrig1);
  TriggerCuts(*(trigPlace[1]), *(trigName[1]), CUTS::eRTrig2);
  updateMet();
  
  // // SET NUMBER OF RECO PARTICLES
  // // MUST BE IN ORDER: Muon/Electron, Tau, Jet
  getGoodRecoLeptons(*_Electron, CUTS::eRElec1, CUTS::eGElec, _Electron->pstats["Elec1"]);
  getGoodRecoLeptons(*_Electron, CUTS::eRElec2, CUTS::eGElec, _Electron->pstats["Elec2"]);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon1, CUTS::eGMuon, _Muon->pstats["Muon1"]);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon2, CUTS::eGMuon, _Muon->pstats["Muon2"]);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau1, CUTS::eGTau, _Tau->pstats["Tau1"]);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau2, CUTS::eGTau, _Tau->pstats["Tau2"]);

  getGoodRecoJets(CUTS::eRJet1, _Jet->pstats["Jet1"]);
  getGoodRecoJets(CUTS::eRJet2, _Jet->pstats["Jet2"]);
  getGoodRecoJets(CUTS::eRCenJet, _Jet->pstats["CentralJet"]);
  getGoodRecoJets(CUTS::eRBJet, _Jet->pstats["BJet"]);
  getGoodRecoJets(CUTS::eR1stJet, _Jet->pstats["FirstLeadingJet"]);
  getGoodRecoJets(CUTS::eR2ndJet, _Jet->pstats["SecondLeadingJet"]);
  ///VBF Susy cut on leadin jets
  VBFTopologyCut(distats["VBFSUSY"]);

  /////lepton lepton topology cuts
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1,CUTS::eRTau1, CUTS::eElec1Tau1, distats["Electron1Tau1"]);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau1, CUTS::eElec2Tau1, distats["Electron2Tau1"]);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1, CUTS::eRTau2, CUTS::eElec1Tau2, distats["Electron1Tau2"]);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau2, CUTS::eElec2Tau2, distats["Electron2Tau2"]);

  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau1, CUTS::eMuon1Tau1, distats["Muon1Tau1"]);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau2, CUTS::eMuon1Tau2, distats["Muon1Tau2"]);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau1, CUTS::eMuon2Tau1, distats["Muon2Tau1"]);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau2, CUTS::eMuon2Tau2, distats["Muon2Tau2"]);

  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec1, CUTS::eMuon1Elec1, distats["Muon1Electron1"]);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec2, CUTS::eMuon1Elec2, distats["Muon1Electron2"]);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec1, CUTS::eMuon2Elec1, distats["Muon2Electron1"]);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec2, CUTS::eMuon2Elec2, distats["Muon2Electron2"]);

  ////DIlepton topology cuts
  getGoodLeptonCombos(*_Tau, *_Tau, CUTS::eRTau1, CUTS::eRTau2, CUTS::eDiTau, distats["DiTau"]);
  getGoodLeptonCombos(*_Electron, *_Electron, CUTS::eRElec1, CUTS::eRElec2, CUTS::eDiElec, distats["DiElectron"]);
  getGoodLeptonCombos(*_Muon, *_Muon, CUTS::eRMuon1, CUTS::eRMuon2, CUTS::eDiMuon, distats["DiMuon"]);

  ////Dijet cuts
  getGoodDiJets(distats["DiJet"]);

  if(event % 50000 == 0) {
  cout << "Event #" << event << endl;
  }
}


////Reads cuts from Cuts.in file and see if the event has enough particles
void Analyzer::fillCuts() {
  const unordered_map<string,pair<int,int> >* cut_info = histo.get_cuts();
  const vector<string>* cut_order = histo.get_cutorder();

  bool prevTrue = true;
  maxCut=0;

  
  for(size_t i = 0; i < cut_order->size(); i++) {
    //  for(vector<string>::const_iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
  
  string cut = cut_order->at(i);
    if(isData && cut.find("Gen") != string::npos) continue;

    int min= cut_info->at(cut).first;
    int max= cut_info->at(cut).second;
    int nparticles = goodParts[cut_num.at(cut)]->size();
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      if((cut_num.at(cut) == CUTS::eR1stJet || cut_num.at(cut) == CUTS::eR2ndJet) && goodParts[cut_num.at(cut)]->at(0) == -1 ) {
	prevTrue = false;
	continue;  ////dirty dirty hack
      }
      cuts_per[i]++;
      cuts_cumul[i] += (prevTrue) ? 1 : 0;
      maxCut += (prevTrue) ? 1 : 0;
    } else prevTrue = false;
  }


}

/////// maxcut made -1 if doesn't pass all of the cuts
void Analyzer::CRfillCuts() {
  const unordered_map<string,pair<int,int> >* cut_info = histo.get_cuts();
  const vector<string>* cut_order = histo.get_cutorder();

  maxCut=0;

  for(size_t i = 0; i < cut_order->size(); i++) {
    string cut = cut_order->at(i);
    if(isData && cut.find("Gen") != string::npos) continue;

    int min= cut_info->at(cut).first;
    int max= cut_info->at(cut).second;
    int nparticles = goodParts[cut_num.at(cut)]->size();
    if( (nparticles < min) || (nparticles > max && max != -1)) {
      maxCut = -1;
      return;
    } else if((cut_num.at(cut) == CUTS::eR1stJet || cut_num.at(cut) == CUTS::eR2ndJet) && goodParts[cut_num.at(cut)]->at(0) == -1 ) {
      maxCut = -1;
      return;
    }
  }
  
  int factor = crbins;
  for(auto tester: testVec) {
    factor /= 2;
    /////get variable value from maper.first.
    if(tester->test(this)) { ///pass cut
      maxCut += factor;
    }
  }
  if(isData && blinded && maxCut == SignalRegion) return;
  cuts_per[maxCut]++;
}


///Prints the number of events that passed each cut per event and cumulatively
//done at the end of the analysis
void Analyzer::printCuts() {
  vector<string> cut_order;
  if(crbins > 1) cut_order = *(histo.get_folders());
  else cut_order = *(histo.get_cutorder());

  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(3);
  cout << "\n";
  cout << "Selection Efficiency " << "\n";
  cout << "Total events: " << nentries << "\n";
  cout << "                        Name                  Indiv.";
  if(crbins == 1) cout << "            Cumulative";
  cout << endl << "---------------------------------------------------------------------------\n";
  for(size_t i = 0; i < cut_order.size(); i++) {
    cout << setw(28) << cut_order.at(i) << "    ";
    if(isData && cut_order.at(i).find("Gen") != string::npos) cout << "Skipped" << endl;
    else if(crbins != 1 && blinded && i == SignalRegion) cout << "Blinded Signal Region" << endl;
    else {
      cout << setw(10) << cuts_per.at(i) << "  ( " << setw(5) << ((float)cuts_per.at(i)) / nentries << ") ";
      if(crbins == 1) cout << setw(12) << cuts_cumul.at(i) << "  ( " << setw(5) << ((float)cuts_cumul.at(i)) / nentries << ") ";

      cout << endl;
    }
  }
  cout << "---------------------------------------------------------------------------\n";  
  histo.fill_histogram();
}

/////////////PRIVATE FUNCTIONS////////////////



///Calculates met from values from each file plus smearing and treating muons as neutrinos
void Analyzer::updateMet() {
  ////// Neutrino update before calculation
  // if(distats["Run"].bmap.at("TreatMuonsAsNeutrinos")) {
  //   for(vec_iter it=goodParts[CUTS::eRMuon1]->begin(); it!=goodParts[CUTS::eRMuon1]->end(); it++) {
  //     if(find(goodParts[CUTS::eRMuon2]->begin(), goodParts[CUTS::eRMuon2]->end(), (*it)) != goodParts[CUTS::eRMuon2]->end() ) continue;
  //     deltaMEx += _Muon->smearP.at(*it).Px();
  //     deltaMEy += _Muon->smearP.at(*it).Py();
  //   }    
  //   for(vec_iter it=goodParts[CUTS::eRMuon2]->begin(); it!=goodParts[CUTS::eRMuon2]->end(); it++) {
  //     deltaMEx += _Muon->smearP.at(*it).Px();
  //     deltaMEy += _Muon->smearP.at(*it).Py();
  //   }
  // }
  ///---MHT and HT calculations----////
  int i=0;
  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it!=_Jet->smearP.end(); it++, i++) {
    if( (it->Pt() > distats["Run"].dmap.at("JetPtForMhtAndHt")) && (fabs(it->Eta()) < distats["Run"].dmap.at("JetEtaForMhtAndHt")) ) {
      if(distats["Run"].bmap.at("ApplyJetLooseIDforMhtAndHt") && !passedLooseJetID(i) ) continue;
      
      sumpxForMht -= it->Px();
      sumpyForMht -= it->Py();
      sumptForHt  += it->Pt();
    }
  }
  phiForMht = atan2(sumpyForMht,sumpxForMht);
  theMETVector.SetPxPyPzE(theMETVector.Px()+deltaMEx, theMETVector.Py()+deltaMEy, theMETVector.Pz(),
  			  TMath::Sqrt(pow(theMETVector.Px()+deltaMEx,2) + pow(theMETVector.Py()+deltaMEy,2)));

  /////MET CUTS
  
  if(!passCutRange("Met", theMETVector.Pt(), distats["Run"])) return;
  if(distats["Run"].bmap.at("DiscrByHT") && sumptForHt < distats["Run"].dmap.at("HtCut")) return; 
  
  goodParts[CUTS::eMET]->push_back(1);
}


/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral(TTree* BOOM, string infile) {
  SetBranch("Trigger_decision", Trigger_decision);
  SetBranch("Trigger_names", Trigger_names);
  SetBranch("nTruePUInteractions", nTruePU);
  SetBranch("bestVertices", bestVertices);
  SetBranch("weightevt", gen_weight);
  SetBranch("Met_type1PF_px", Met[0]);
  SetBranch("Met_type1PF_py", Met[1]);
  SetBranch("Met_type1PF_pz", Met[2]);

  read_info(FILESPACE + "ElectronTau_info.in");
  read_info(FILESPACE + "MuonTau_info.in");
  read_info(FILESPACE + "MuonElectron_info.in");
  read_info(FILESPACE + "DiParticle_info.in");
  read_info(FILESPACE + "VBFCuts_info.in");
  read_info(FILESPACE + "Run_info.in");

  BOOM->GetEntry(0);
  for(int i = 0; i < nTrigReq; i++) {
    for(int j = 0; j < (int)trigName[i]->size(); j++) {
      for(int k = 0; k < (int)Trigger_names->size(); k++) {
	if(Trigger_names->at(k).find(trigName[i]->at(j)) != string::npos) {
	  trigPlace[i]->at(j) = k;
	  break;
	}
      }
    }
  }
  BOOM->SetBranchStatus("Trigger_names", 0);

}


///parsing method that gets info on diparts and basic run info
//put in map called "distats"
void Analyzer::read_info(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }


  string group, line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    vector<string> stemp;

    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      cout << "error in " << filename << "; no groups specified for data" << endl;
      exit(1);
    } else if(stemp.size() == 2) {
      if(stemp.at(0).find("Trigger") != string::npos) {
	int ntrig = (stemp.at(0).find("1") != string::npos) ? 0 : 1;
	trigName[ntrig]->push_back(stemp.at(1));
	trigPlace[ntrig]->push_back(0);
	continue;
      }
	
      char* p;
      strtod(stemp[1].c_str(), &p);
      if(stemp[1] == "1" || stemp[1] == "true") distats[group].bmap[stemp[0]]=true;
      else if(stemp[1] == "0" || stemp[1] == "false") distats[group].bmap[stemp[0]]=false; 
      else if(*p) distats[group].smap[stemp[0]] = stemp[1];
      else  distats[group].dmap[stemp[0]]=stod(stemp[1]);

    } else if(stemp.size() == 3) distats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
  }
  info_file.close();
}


// This code works pretty much (at least in my tests), but dagnabit, its ugly.  They all can't be winners, at least now...
void Analyzer::setCutNeeds() {

  for(auto e: Enum<CUTS>()) {
    need_cut[e] = false;
  }
  for(auto it: *histo.get_groups()) {
    if(fillInfo[it]->type == FILLER::None) continue;
    need_cut[fillInfo[it]->ePos] = true;
    if(adjList.find(fillInfo[it]->ePos) == adjList.end()) continue;
    for(auto e: adjList.at(fillInfo[it]->ePos)) {
      need_cut[e] = true;
    }
  }
  for(auto it: *histo.get_cutorder()) {
    need_cut[cut_num.at(it)] = true;
    if(adjList.find(cut_num.at(it)) == adjList.end()) continue;
    for(auto e: adjList.at(cut_num.at(it))) {
      need_cut[e] = true;
    }
  }
  for(auto it: testVec) {
    CUTS ePos = it->info->ePos;
    need_cut[ePos] = true;
    if(adjList.find(ePos) == adjList.end()) continue;
    for(auto e: adjList.at(ePos)) {
      need_cut[e] = true;
    }
  }

  for(auto it: _Jet->extraCuts) {
    need_cut[it] = true;
    if(adjList.find(it) == adjList.end()) continue;
    for(auto e: adjList.at(it)) {
      need_cut[e] = true;
    }
  }


  for(auto it: jetCuts) {
    if(need_cut[it]) {
      for(auto it2: _Jet->overlapCuts(it)) {
	need_cut[it2] = true;
	if(adjList.find(it2) == adjList.end()) continue;
	for(auto e: adjList.at(it2)) {
	  need_cut[e] = true;
	}
      }
    }
  }


  if( !(need_cut[CUTS::eRTau1] || need_cut[CUTS::eRTau2]) ) {
    _Tau->unBranch();
  } else {
    for(auto it: _Tau->extraCuts) {
      need_cut[it] = true;
      if(adjList.find(it) == adjList.end()) continue;
      for(auto e: adjList.at(it)) {
	need_cut[e] = true;
      }
    }
  }
  if( !(need_cut[CUTS::eRElec1] || need_cut[CUTS::eRElec2]) ) {
    _Electron->unBranch();
  } else {
    for(auto it: _Electron->extraCuts) {
      need_cut[it] = true;
      if(adjList.find(it) == adjList.end()) continue;
      for(auto e: adjList.at(it)) {
	need_cut[e] = true;
      }
    }
  }
  if( !(need_cut[CUTS::eRMuon1] || need_cut[CUTS::eRMuon2]) ) {
    _Muon->unBranch();
  } else {
    for(auto it: _Muon->extraCuts) {
      need_cut[it] = true;
      if(adjList.find(it) == adjList.end()) continue;
      for(auto e: adjList.at(it)) {
	need_cut[e] = true;
      }
    }
  }

  if(isData) return;

  bool passGen = false;
  for(auto e: genCuts) {
    passGen = passGen || need_cut[e];
   
  }
  if(!passGen) _Gen->unBranch();
  else {
    if(need_cut[CUTS::eGTau]) genMaper[15] = new GenFill(2, CUTS::eGTau);
    if(need_cut[CUTS::eGTop]) genMaper[6] = new GenFill(2, CUTS::eGTop);
    if(need_cut[CUTS::eGElec]) genMaper[11] = new GenFill(1, CUTS::eGElec);
    if(need_cut[CUTS::eGMuon]) genMaper[13] = new GenFill(1, CUTS::eGMuon);
    if(need_cut[CUTS::eGZ]) genMaper[23] = new GenFill(2, CUTS::eGZ);
    if(need_cut[CUTS::eGW]) genMaper[24] = new GenFill(2, CUTS::eGW);
    if(need_cut[CUTS::eGHiggs]) genMaper[25] = new GenFill(2, CUTS::eGHiggs);
    //  , CUTS::eNuTau
  }
  
}


///Smears lepton only if specified and not a data file.  Otherwise, just filles up lorentz vectors
//of the data into the vector container smearP with is in each lepton object.
void Analyzer::smearLepton(Lepton& lepton, CUTS eGenPos, const PartStats& stats) {
  lepton.smearP.clear();

  for(int i = 0; i < (int)lepton.pt->size(); i++) {
    TLorentzVector tmpSmear;
    tmpSmear.SetPtEtaPhiE(lepton.pt->at(i), lepton.eta->at(i), lepton.phi->at(i), lepton.energy->at(i));

    if(isData || !stats.bmap.at("SmearTheParticle")) {
      lepton.smearP.push_back(tmpSmear);
      continue;
    }

    TLorentzVector genVec =  matchLeptonToGen(tmpSmear, lepton.pstats["Smear"],eGenPos);
    if(genVec == TLorentzVector(0,0,0,0)) {      
      lepton.smearP.push_back(tmpSmear);
      continue;
    }

    double smearedPt = (genVec.Pt()*stats.dmap.at("PtScaleOffset")) + (tmpSmear.Pt() - genVec.Pt())*stats.dmap.at("PtSigmaOffset");
    double smearedEta =(genVec.Eta()*stats.dmap.at("EtaScaleOffset")) + (tmpSmear.Eta() - genVec.Eta())*stats.dmap.at("EtaSigmaOffset");
    double smearedPhi = (genVec.Phi() * stats.dmap.at("PhiScaleOffset")) + (tmpSmear.Phi() - genVec.Phi())*stats.dmap.at("PhiSigmaOffset");
    double smearedEnergy = (genVec.Energy()*stats.dmap.at("EnergyScaleOffset")) + (tmpSmear.Energy() - genVec.Energy())*stats.dmap.at("EnergySigmaOffset");
    
    TLorentzVector final;
    final.SetPtEtaPhiE(smearedPt, smearedEta, smearedPhi, smearedEnergy);
    lepton.smearP.push_back(final);
    deltaMEx += tmpSmear.Px() - final.Px();
    deltaMEy += tmpSmear.Py() - final.Py();

  }
 
}

///Same as smearlepton, just jet specific
void Analyzer::smearJet(const PartStats& stats) {
  _Jet->smearP.clear();

  TLorentzVector jetV;

  for(size_t i=0; i< _Jet->pt->size(); i++) {
    jetV.SetPtEtaPhiE(_Jet->pt->at(i), _Jet->eta->at(i), _Jet->phi->at(i), _Jet->energy->at(i));

    if(isData || !stats.bmap.at("SmearTheJet")) {
      _Jet->smearP.push_back(jetV);
      continue;
    }
    
    if(JetMatchesLepton(*_Muon, jetV, stats.dmap.at("MuonMatchingDeltaR"), CUTS::eGMuon) ||
       JetMatchesLepton(*_Tau, jetV, stats.dmap.at("TauMatchingDeltaR"), CUTS::eGTau) ||       
       JetMatchesLepton(*_Electron, jetV,stats.dmap.at("ElectronMatchingDeltaR"), CUTS::eGElec)){

      _Jet->smearP.push_back(jetV);

      continue;
    }

    _Jet->smearP.push_back(stats.dmap.at("JetEnergyScaleOffset") * jetV);
    deltaMEx += (1 - stats.dmap.at("JetEnergyScaleOffset"))*jetV.Px();
    deltaMEy += (1 -stats.dmap.at("JetEnergyScaleOffset"))*jetV.Py();
  }
}

/////checks if jet is close to a lepton and the lepton is a gen particle, then the jet is a lepton object, so
//this jet isn't smeared
bool Analyzer::JetMatchesLepton(const Lepton& lepton, const TLorentzVector& jetV, double partDeltaR, CUTS eGenPos) {
  TLorentzVector tempV;
  for(size_t j = 0; j < lepton.pt->size(); j++) {
    tempV.SetPtEtaPhiE(lepton.pt->at(j), lepton.eta->at(j), lepton.phi->at(j), lepton.energy->at(j));
    if(jetV.DeltaR(tempV) < partDeltaR && matchLeptonToGen(tempV, lepton.pstats.at("Smear"), eGenPos) != TLorentzVector(0,0,0,0)) return true;
  }
  return false;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchLeptonToGen(const TLorentzVector& lvec, const PartStats& stats, CUTS ePos) {
  if(ePos == CUTS::eGTau) {
    return matchTauToGen(lvec, stats.dmap.at("GenMatchingDeltaR"));
  }
  TLorentzVector genVec = TLorentzVector(0,0,0,0);
  
  for(vec_iter it=goodParts[ePos]->begin(); it !=goodParts[ePos]->end();it++) {
    genVec.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
    if(lvec.DeltaR(genVec) <= stats.dmap.at("GenMatchingDeltaR")) {
      if(stats.bmap.at("UseMotherID") && abs(_Gen->motherpdg_id->at(*it)) != stats.dmap.at("MotherID")) continue; 
      return genVec;
    }
  }
  
  return TLorentzVector(0,0,0,0);
}


///Tau specific matching fucntion.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically
TLorentzVector Analyzer::matchTauToGen(const TLorentzVector& lvec, double lDeltaR) {
  TLorentzVector genVec(0,0,0,0);
  int i = 0;
  for(vec_iter it=goodParts[CUTS::eGTau]->begin(); it !=goodParts[CUTS::eGTau]->end();it++, i++) {
    int nu = goodParts[CUTS::eNuTau]->at(i);
    if(nu == -1) continue;

    TLorentzVector tmp1, tmp2;
    tmp1.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
    tmp2.SetPtEtaPhiE(_Gen->pt->at(nu), _Gen->eta->at(nu), _Gen->phi->at(nu), _Gen->energy->at(nu));
    genVec = tmp1 - tmp2;
    if(lvec.DeltaR(genVec) <= lDeltaR) {
      return genVec;
    }
  }
  return TLorentzVector(0,0,0,0);

}


////Calculates the number of gen particles.  Based on id number and status of each particle
void Analyzer::getGoodGen(const PartStats& stats) {

  for(size_t j = 0; j < _Gen->pt->size(); j++) {
    int id = abs(_Gen->pdg_id->at(j)); 
    if(genMaper[id] != nullptr && _Gen->pdg_id->at(j) == genMaper[id]->status) {
      if(id == 15 && (_Gen->pt->at(j) < stats.pmap.at("TauPtCut").first || _Gen->pt->at(j) > stats.pmap.at("TauPtCut").second || abs(_Gen->eta->at(j)) > stats.dmap.at("TauEtaCut"))) continue;
      goodParts[genMaper[id]->ePos]->push_back(j);
    }
  }

    
    // if((abs(_Gen->pdg_id->at(j)) == particle_id) && (_Gen->status->at(j) == particle_status)) {
    //   goodParts[ePos]->push_back(j);
    // }
  //}

}

////Tau neutrino specific function used for calculating the number of hadronic taus
void Analyzer::getGoodTauNu() {

  for(vec_iter it=goodParts[CUTS::eGTau]->begin(); it !=goodParts[CUTS::eGTau]->end();it++) {
    bool leptonDecay = false;
    int nu = -1;
    for(size_t j = 0; j < _Gen->pt->size(); j++) {
      if(abs(_Gen->BmotherIndex->at(j)) == (*it)) {
	if( (abs(_Gen->pdg_id->at(j)) == 16) && (abs(_Gen->motherpdg_id->at(j)) == 15) && (_Gen->status->at(_Gen->BmotherIndex->at(j)) == 2) ) nu = j;
	else if( (abs(_Gen->pdg_id->at(j)) == 12) || (abs(_Gen->pdg_id->at(j)) == 14) ) leptonDecay = true;
      }
    }
    nu = (leptonDecay) ? -1 : nu;
    goodParts[CUTS::eNuTau]->push_back(nu);
  }

}

///Function used to find the number of reco leptons that pass the various cuts.
///Divided into if blocks for the different lepton requirements.
void Analyzer::getGoodRecoLeptons(const Lepton& lep, const CUTS ePos, const CUTS eGenPos, const PartStats& stats) {
  if(! need_cut[ePos]) return;
  int i = 0;

  for(vector<TLorentzVector>::const_iterator it=lep.smearP.begin(); it != lep.smearP.end(); it++, i++) {
    TLorentzVector lvec = (*it);

    if (fabs(lvec.Eta()) > stats.dmap.at("EtaCut")) continue;
    if (lvec.Pt() < stats.pmap.at("PtCut").first || lvec.Pt() > stats.pmap.at("PtCut").second) continue;

    if((lep.pstats.at("Smear").bmap.at("MatchToGen")) && (!isData)) {   /////check
      if(matchLeptonToGen(lvec, lep.pstats.at("Smear") ,eGenPos) == TLorentzVector(0,0,0,0)) continue;
    }

    if (stats.bmap.at("DoDiscrByIsolation")) {
      double firstIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").first : ival(ePos) - ival(CUTS::eRTau1) + 1;
      double secondIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").second : 0;
      if(!lep.get_Iso(i, firstIso, secondIso)) continue;
    }

    if ((lep.type != PType::Tau) && stats.bmap.at("DiscrIfIsZdecay")) { 
      if(isZdecay(lvec, lep)) continue;
    }
    if(!passCutRange("MetDphi", absnormPhi(lvec.Phi() - theMETVector.Phi()), stats)) continue;
    if(!passCutRange("MetMt", calculateLeptonMetMt(lvec), stats)) continue;


    if(lep.type == PType::Muon) {      ////////////////MUON CUTS/////////////
      if(stats.bmap.at("DoDiscrByTightID") && (_Muon->tight->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrBySoftID") && (_Muon->soft->at(i) == 0)) continue;
      
    } else if(lep.type == PType::Electron) {    ///////////////ELECTRON CUT///////////

      //----Require electron to pass ID discriminators
      if(stats.bmap.at("DoDiscrByVetoID") && (_Electron->isPassVeto->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByLooseID") && (_Electron->isPassLoose->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByMediumID") && (_Electron->isPassMedium->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByTightID") && (_Electron->isPassTight->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByHEEPID") && (_Electron->isPassHEEPId->at(i) == 0)) continue;

     
    } else if(lep.type == PType::Tau) {   /////////////TAU CUT/////////////////
      if (stats.bmap.at("DoDiscrByLeadTrack")) {
	if(_Tau->leadChargedCandPt->at(i) < stats.dmap.at("LeadTrackThreshold")) continue;
      }
     
      // ----Require 1 or 3 prongs
      if(stats.smap.at("DiscrByProngType").find("hps") != string::npos && _Tau->decayModeFindingNewDMs->at(i) == 0) continue;
      if(!passProng(stats.smap.at("DiscrByProngType"), _Tau->nProngs->at(i))) continue;

      // ----Electron and Muon vetos
      vector<int>* against = (ePos == CUTS::eRTau1) ? _Tau->againstElectron.first : _Tau->againstElectron.second;
      if (stats.bmap.at("DoDiscrAgainstElectron") && against->at(i) == 0) continue;
      else if (stats.bmap.at("SelectTausThatAreElectrons") && against->at(i) > 0) continue;

      against = (ePos == CUTS::eRTau1) ? _Tau->againstMuon.first : _Tau->againstMuon.second;
      if (stats.bmap.at("DoDiscrAgainstMuon") && against->at(i) == 0) continue;
      else if (stats.bmap.at("SelectTausThatAreMuons") && against->at(i) > 0) continue;

      if (stats.bmap.at("DoDiscrByCrackCut") && isInTheCracks(lvec.Eta())) continue;

      // ----anti-overlap requirements
      if (stats.bmap.at("RemoveOverlapWithMuon1s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithMuon2s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithElectron1s") && isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithElectron2s") && isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"))) continue;
    }
    goodParts[ePos]->push_back(i);    
  }
  
}

////Jet specific function for finding the number of jets that pass the cuts.
//used to find the nubmer of good jet1, jet2, central jet, 1st and 2nd leading jets and bjet.
void Analyzer::getGoodRecoJets(CUTS ePos, const PartStats& stats) {
  if(! need_cut[ePos]) return;
  int i=0;
  

  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it != _Jet->smearP.end(); it++, i++) {
    TLorentzVector lvec = (*it);
    ///if else loop for central jet requirements

    if( ePos == CUTS::eRCenJet) {
      if(fabs(lvec.Eta()) > 2.5) continue;
    } else if (fabs(lvec.Eta()) < stats.pmap.at("EtaCut").first || fabs(lvec.Eta()) > stats.pmap.at("EtaCut").second) continue;

    if (lvec.Pt() < stats.dmap.at("PtCut")) continue;

    /// BJet specific
    if(ePos == CUTS::eRBJet) {
      if(stats.bmap.at("ApplyJetBTagging") && _Jet->bDiscriminator->at(i) <= stats.dmap.at("JetBTaggingCut")) continue;
      if((stats.bmap.at("MatchBToGen")) && !isData && abs(_Jet->partonFlavour->at(i)) != 5) continue;
    } else if (stats.bmap.at("ApplyLooseID") && !passedLooseJetID(i)) continue;

  // ----anti-overlap requirements
    if(stats.bmap.at("RemoveOverlapWithMuon1s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithMuon2s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithElectron1s") && isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithElectron2s") && isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithTau1s") && isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithTau2s") && isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"))) continue;

    /////fill up array
    goodParts[ePos]->push_back(i);    
  }

  if(ePos == CUTS::eR1stJet || ePos == CUTS::eR2ndJet) {
    int potential = -1;
    double prevPt = -1;
    for(vec_iter leadit = goodParts[ePos]->begin(); leadit != goodParts[ePos]->end(); ++leadit) {
      if(((ePos == CUTS::eR2ndJet && (*leadit) != leadIndex) || ePos == CUTS::eR1stJet) && _Jet->smearP.at(*leadit).Pt() > prevPt) {
	potential = (*leadit);
	prevPt = _Jet->smearP.at(*leadit).Pt();
      }
    }
    goodParts[ePos]->clear();
    goodParts[ePos]->push_back(potential);
    if(ePos == CUTS::eR1stJet) leadIndex = goodParts[CUTS::eR1stJet]->at(0); 
  }

} 

///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(vec_iter it=goodParts[ePos]->begin(); it != goodParts[ePos]->end(); it++) {
    if(lvec.DeltaR(overlapper.smearP.at(*it)) < MatchingDeltaR) return true;
  }
  return false;
}

///Tests if tau decays into the specified number of jet prongs.
bool Analyzer::passProng(string prong, int value) {
  return ( (prong.find("1") != string::npos && value == 1) ||
	   (prong.find("2") != string::npos && value == 2) ||
	   (prong.find("3") != string::npos && value == 3) );
}


////Tests if tau is within the cracks of the detector (the specified eta ranges)
bool Analyzer::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 ||
          (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
          (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
          (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
          (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}
 

//Tests if a jet meets a litany of different tests
bool Analyzer::passedLooseJetID(int nobj) {
  if (_Jet->neutralHadEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->neutralEmEmEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->numberOfConstituents->at(nobj) <= 1) return false;
  if (_Jet->muonEnergyFraction->at(nobj) >= 0.80) return false;
  if ( (fabs(_Jet->smearP.at(nobj).Eta()) < 2.4) && 
           ((_Jet->chargedHadronEnergyFraction->at(nobj) <= 0.0) || 
	    (_Jet->chargedMultiplicity->at(nobj) <= 0.0) || 
	    (_Jet->chargedEmEnergyFraction->at(nobj) >= 0.99) )) return false;
  return true;
}


///sees if the event passed one of the two cuts provided
void Analyzer::TriggerCuts(vector<int>& prevTrig, const vector<string>& trigvec, CUTS ePos) {
  if(! need_cut[ePos]) return;
  for(int i = 0; i < (int)trigvec.size(); i++) {
    if(Trigger_decision->at(prevTrig.at(i)) == 1) {
      goodParts[ePos]->push_back(0);
      return;
    }
  }
}


////VBF specific cuts dealing with the leading jets.
void Analyzer::VBFTopologyCut(const PartStats& stats) {
  if(! need_cut[CUTS::eSusyCom]) return;

  if(goodParts[CUTS::eR1stJet]->at(0) == -1 || goodParts[CUTS::eR2ndJet]->at(0) == -1) return;

  TLorentzVector ljet1 = _Jet->smearP.at(goodParts[CUTS::eR1stJet]->at(0));
  TLorentzVector ljet2 = _Jet->smearP.at(goodParts[CUTS::eR2ndJet]->at(0));
  
  if(!passCutRange("Mass", (ljet1 + ljet2).M(), stats)) return;
  if(!passCutRange("Pt", (ljet1 + ljet2).Pt(), stats)) return;
  if(!passCutRange("DeltaEta", abs(ljet1.Eta() - ljet2.Eta()), stats)) return;
  if(!passCutRange("DeltaEta", absnormPhi(ljet1.Phi() - ljet2.Phi()), stats)) return;
  
  if(stats.bmap.at("DiscrByOSEta")) {
    if((ljet1.Eta() * ljet2.Eta()) >= 0) return;
  }

  double dphi1 = normPhi(ljet1.Phi() - theMETVector.Phi());
  double dphi2 = normPhi(ljet2.Phi() - theMETVector.Phi());
  double r1, r2, alpha;
  
  if(stats.bmap.at("DiscrByR1")) {
    r1 = sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) );
    if(r1 < stats.pmap.at("R1Cut").first || r1 > stats.pmap.at("R1Cut").second) return;

  }
  if(stats.bmap.at("DiscrByR2")) {
    r2 = sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0) );
    if(r2 < stats.pmap.at("R2Cut").first || r2 > stats.pmap.at("R2Cut").second) return;
  }
  if(stats.bmap.at("DiscrByAlpha")) {
    TLorentzVector addVec = ljet1 + ljet2;
    alpha = (addVec.M() > 0) ? ljet2.Pt() / addVec.M() : -1;
    if(alpha < stats.pmap.at("AlphaCut").first || alpha > stats.pmap.at("AlphaCut").second) return;
  }
  if( !passCutRange("Dphi1", abs(dphi1), stats)) return;
  if( !passCutRange("Dphi2", abs(dphi2), stats)) return;

  goodParts[CUTS::eSusyCom]->push_back(0);
}

inline bool Analyzer::passCutRange(string CutName, double value, const PartStats& stats) {
  return ( !(stats.bmap.at("DiscrBy" + CutName)) || (value > stats.pmap.at(CutName + "Cut").first && value < stats.pmap.at(CutName + "Cut").second) );

}


//-----Calculate lepton+met transverse mass
double Analyzer::calculateLeptonMetMt(const TLorentzVector& Tobj) {
  double px = Tobj.Px() + theMETVector.Px();
  double py = Tobj.Py() + theMETVector.Py();
  double et = Tobj.Et() + theMETVector.Energy(); //TMath::Sqrt((theMETVector.Px() * theMETVector.Px()) + (theMETVector.Py() * theMETVector.Py()));
  double mt2 = et*et - (px*px + py*py);
  return (mt2 >= 0) ? sqrt(mt2) : -1;
}


/////Calculate the diparticle mass based on how to calculate it
///can use Collinear Approximation, which can fail (number failed available in a histogram)
///can use VectorSumOfVisProductAndMet which is sum of particles and met
///Other which is adding without met
double Analyzer::diParticleMass(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string howCalc) {
  bool ratioNotInRange = false;
  TLorentzVector The_LorentzVect;


  //////check this equation/////
  if(howCalc == "CollinearApprox") {
    double denominator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1 = (Tobj2.Py()*theMETVector.Px() - Tobj2.Px()*theMETVector.Py())/denominator;
    double x2 = (Tobj1.Px()*theMETVector.Py() - Tobj1.Py()*theMETVector.Px())/denominator;
    ratioNotInRange=!((x1 < 0.) && (x2 < 0.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px()*(1 + x1) + Tobj2.Px()*(1+x2)), (Tobj1.Py()*(1+x1) + Tobj2.Py()*(1+x2)), (Tobj1.Pz()*(1+x1) + Tobj2.Pz()*(1+x2)), (Tobj1.Energy()*(1+x1) + Tobj2.Energy()*(1+x2)) );
      return The_LorentzVect.M();
    }
  } 

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    return (Tobj1 + Tobj2 + theMETVector).M();
  }

  return (Tobj1 + Tobj2).M();
}

////Tests if the CollinearApproximation works for finding the mass of teh particles
bool Analyzer::passDiParticleApprox(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string howCalc) {
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + theMETVector.Px())) - (Tobj2.Px() * (Tobj1.Py() + theMETVector.Py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + theMETVector.Py())) - (Tobj1.Py() * (Tobj2.Px() + theMETVector.Px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    return (x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.);
  } else {
    return true;
  }
}




/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonCombos(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats) {
  if(! need_cut[ePosFin]) return;
  bool sameParticle = (&lep1 == &lep2);
  TLorentzVector part1, part2;

  for(vec_iter i1=goodParts[ePos1]->begin(); i1 != goodParts[ePos1]->end(); i1++) {
    for(vec_iter i2=goodParts[ePos2]->begin(); i2 != goodParts[ePos2]->end(); i2++) {
      if(sameParticle && (*i2) <= (*i1)) continue;
      part1 = lep1.smearP.at(*i1);
      part2 = lep2.smearP.at(*i2);

      if(stats.bmap.at("DiscrByDeltaR") && (part1.DeltaR(part2)) < stats.dmap.at("DeltaRCut")) continue;
      if(stats.smap.at("DiscrByOSLSType") == "LS" && (lep1.charge->at(*i1) * lep2.charge->at(*i2) <= 0)) continue;
      else if(stats.smap.at("DiscrByOSLSType") == "OS" && (lep1.charge->at(*i1) * lep2.charge->at(*i2) >= 0)) continue;

      if( !passCutRange("CosDphi", cos(absnormPhi( part1.Phi() - part2.Phi())), stats)) continue;

  // ----Mass window requirement
      
      if (stats.bmap.at("DiscrByMassReco")) {
      	double diMass = diParticleMass(part1,part2, stats.smap.at("HowCalculateMassReco"));
      	if( diMass < stats.pmap.at("MassCut").first || diMass > stats.pmap.at("MassCut").second) continue;
      }

      if (stats.bmap.at("DiscrByCDFzeta2D")) {
      	double CDFzeta = stats.dmap.at("PZetaCutCoefficient") * getPZeta(part1, part2).first 
	  + stats.dmap.at("PZetaVisCutCoefficient") * getPZeta(part1, part2).second;
      	if( CDFzeta < stats.pmap.at("CDFzeta2DCutValue").first || CDFzeta > stats.pmap.at("CDFzeta2DCutValue").second ) continue;
      }

      ///      optional cut
      // if(stats.bmap.find("DeltaPtAndMet") != stats.bmap.end() && stats.bmap.at("DiscrByCosDphi_DeltaPtAndMet")) {
      // 	double DPhi = absnormPhi(atan2(part1.Py() + part2.Py(), part1.Px() + part2.Px()) - theMETVector.Phi());
      // 	if( cos(DPhi) < stats.pmap.at("CosDphi_DeltaPtMetCut").first || cos(DPhi) > stats.pmap.at("CosDphi_DeltaPtMetCut").second) continue;
      // }

      //////////abs on the difference????
      ///////////////////

      if (stats.bmap.at("DiscrByDeltaPtDivSumPt")) {
	double ptDiv = (lep1.smearP.at(*i1).Pt() - lep2.smearP.at(*i2).Pt()) / (lep1.smearP.at(*i1).Pt() + lep2.smearP.at(*i2).Pt());
	if( ptDiv < stats.pmap.at("DeltaPtDivSumPtCutValue").first || ptDiv > stats.pmap.at("DeltaPtDivSumPtCutValue").second) continue;
      }

      if (stats.bmap.at("DiscrByDeltaPt")) {
	double deltaPt = lep1.smearP.at(*i1).Pt() - lep2.smearP.at(*i2).Pt(); 
	if(deltaPt < stats.pmap.at("DeltaPtCutValue").first || deltaPt > stats.pmap.at("DeltaPtCutValue").second) continue;
      }

      ///Particlesp that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2 
      goodParts[ePosFin]->push_back((*i1)*BIG_NUM + (*i2));
    }
  }
}


//////////////LOOK INTO DIJET PICKING
///////HOW TO GET RID OF REDUNCENCIES??

/////Same as gooddilepton, just jet specific
void Analyzer::getGoodDiJets(const PartStats& stats) {
  if(! need_cut[CUTS::eDiJet]) return;
   TLorentzVector jet1, jet2;
  // ----Separation cut between jets (remove overlaps)
  for(vec_iter ij2=goodParts[CUTS::eRJet2]->begin(); ij2 != goodParts[CUTS::eRJet2]->end(); ij2++) {
    jet2 = _Jet->smearP.at(*ij2);
    for(vec_iter ij1=goodParts[CUTS::eRJet1]->begin(); ij1 != goodParts[CUTS::eRJet1]->end() && (*ij1) < (*ij2); ij1++) {
      jet1 = _Jet->smearP.at(*ij1);

      if (stats.bmap.at("DiscrByDeltaR")) {
	if(jet1.DeltaR(jet2) < stats.dmap.at("DeltaRCut")) continue;
      }

      if( !passCutRange("DeltaEta", abs(jet1.Eta() - jet2.Eta()), stats) ) continue;
      if( !passCutRange("DeltaPhi", abs(jet1.Phi() - jet2.Phi()), stats) ) continue;

      if (stats.bmap.at("DiscrByOSEta")) {
	if((jet1.Eta() * jet2.Eta()) >= 0) continue;
      }
      // ----Require both legs to be almost back-to-back in phi
      if( !passCutRange("CosDphi", cos(absnormPhi(jet1.Phi() - jet2.Phi())), stats) ) continue;

      // ----Mass window requirement
      if (stats.bmap.at("DiscrByMassReco")) {
	if( ((jet1 + jet2).M() < stats.pmap.at("MassCut").first) || ((jet1 + jet2).M() > stats.pmap.at("MassCut").second) ) continue;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2 
      goodParts[CUTS::eDiJet]->push_back((*ij1)*_Jet->smearP.size() + (*ij2));
    }
  }
}


///////Only tested for if is Zdecay, can include massptasymmpair later?
/////Tests to see if a light lepton came form a zdecay
bool Analyzer::isZdecay(const TLorentzVector& theObject, const Lepton& lep) {
  bool eventIsZdecay = false;
  const float zMass = 90.1876;
  const float zWidth = 2.4952;
  float zmmPtAsymmetry = -10.;

  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.
  for(vector<TLorentzVector>::const_iterator lepit= lep.smearP.begin(); lepit != lep.smearP.end(); lepit++) {
    if(theObject.DeltaR(*lepit) < 0.3) continue;
    if(theObject == (*lepit)) continue;

    TLorentzVector The_LorentzVect = theObject + (*lepit);
    zmmPtAsymmetry = (theObject.Pt() - lepit->Pt()) / (theObject.Pt() + lepit->Pt());

    if( (abs(The_LorentzVect.M() - zMass) < 3.*zWidth) || (fabs(zmmPtAsymmetry) < 0.20) ) {
      eventIsZdecay = true;
      break;
    }
  }

  return eventIsZdecay;
}


///Calculates the Pzeta value
pair<double, double> Analyzer::getPZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double px = visPx + theMETVector.Px();
  double py = visPy + theMETVector.Py();
  return make_pair(px*zetaX + py*zetaY, visPx*zetaX + visPy*zetaY);
}


///Normalizes phi to be between -PI and PI

///Takes the absolute value of of normPhi (made because constant use)



////Grabs a list of the groups of histograms to be filled and asked Fill_folder to fill up the histograms
void Analyzer::fill_histogram() {
  if(distats["Run"].bmap["ApplyGenWeight"] && gen_weight == 0.0) return;

  if(crbins != 1) CRfillCuts();
  else fillCuts();
  if(isData && blinded && maxCut == SignalRegion) return;
  const vector<string>* groups = histo.get_groups();
  wgt = pu_weight;
  if(distats["Run"].bmap["ApplyGenWeight"]) wgt *= (gen_weight > 0) ? 1.0 : -1.0;

  for(auto it: *groups) {
    fill_Folder(it, maxCut);
  }
}

///Function that fills up the histograms
void Analyzer::fill_Folder(string group, const int max) {
  if(group == "FillRun") {
    if(crbins != 1) {
      for(int i = 0; i < crbins; i++) {
	histo.addVal(false, group, i, "Events", 1);
      }
    }  
    else histo.addVal(false, group,histo.get_maxfolder(), "Events", 1);
    histAddVal(true, "Events");
    histAddVal(bestVertices, "NVertices");
  } else if(!isData && group == "FillGen") {

    int nhadtau = 0;
    TLorentzVector genVec;
    int i = 0;
    for(vec_iter it=goodParts[CUTS::eGTau]->begin(); it!=goodParts[CUTS::eGTau]->end(); it++, i++) {

      int nu = goodParts[CUTS::eNuTau]->at(i);
      if(nu != -1) {
	TLorentzVector tmp1, tmp2;
	tmp1.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
	tmp2.SetPtEtaPhiE(_Gen->pt->at(nu), _Gen->eta->at(nu), _Gen->phi->at(nu), _Gen->energy->at(nu));
	genVec = tmp1 - tmp2;
	histAddVal(genVec.Pt(), "HadTauPt");
	histAddVal(genVec.Eta(), "HadTauEta");
	nhadtau++;

      }
      histAddVal(_Gen->energy->at(*it), "TauEnergy");
      histAddVal(_Gen->pt->at(*it), "TauPt");
      histAddVal(_Gen->eta->at(*it), "TauEta");
      histAddVal(_Gen->phi->at(*it), "TauPhi");
      for(vec_iter it2=it+1; it2!=goodParts[CUTS::eGTau]->end(); it2++) {
	TLorentzVector genObjt1;
	genObjt1.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
	TLorentzVector genObjt2;
	genObjt2.SetPtEtaPhiE(_Gen->pt->at(*it2), _Gen->eta->at(*it2), _Gen->phi->at(*it2), _Gen->energy->at(*it2));
	histAddVal(diParticleMass(genObjt1,genObjt2, "none"), "DiTauMass");
      }
    }
    histAddVal(goodParts[CUTS::eGTau]->size(), "NTau");
    histAddVal(nhadtau, "NHadTau");

    for(vec_iter it=goodParts[CUTS::eGMuon]->begin(); it!=goodParts[CUTS::eGMuon]->end(); it++) {
      histAddVal(_Gen->energy->at(*it), "MuonEnergy");
      histAddVal(_Gen->pt->at(*it), "MuonPt");
      histAddVal(_Gen->eta->at(*it), "MuonEta");
      histAddVal(_Gen->phi->at(*it), "MuonPhi");
    }
    histAddVal(goodParts[CUTS::eGMuon]->size(), "NMuon");

          

  } else if(fillInfo[group]->type == FILLER::Single) {
    Particle* part = fillInfo[group]->part;
    CUTS ePos = fillInfo[group]->ePos;

    for(vec_iter it=goodParts[ePos]->begin(); it!=goodParts[ePos]->end(); it++) {
      histAddVal(part->smearP.at(*it).Energy(), "Energy");
      histAddVal(part->smearP.at(*it).Pt(), "Pt");
      histAddVal(part->smearP.at(*it).Eta(), "Eta");
      histAddVal(part->smearP.at(*it).Phi(), "Phi");
      if(part->type == PType::Tau) {
	histAddVal(_Tau->nProngs->at(*it), "NumSignalTracks");
  	histAddVal(_Tau->charge->at(*it), "Charge");
	histAddVal(_Tau->leadChargedCandPt->at(*it), "SeedTrackPt");
      } 
      if(part->type != PType::Jet) {
      	histAddVal(calculateLeptonMetMt(part->smearP.at(*it)), "MetMt");  
      }
    }
    if((part->type != PType::Jet ) && goodParts[ePos]->size() > 0) {
      double leadpt = 0;
      double leadeta = 0;
      for(vec_iter it=goodParts[ePos]->begin(); it!=goodParts[ePos]->end(); it++) {
    	if(part->smearP.at(*it).Pt() >= leadpt) {
    	  leadpt = part->smearP.at(*it).Pt();
    	  leadeta = part->smearP.at(*it).Eta();
    	}
      }

      histAddVal(leadpt, "FirstLeadingPt");
      histAddVal(leadeta, "FirstLeadingEta");
    }

    histAddVal(goodParts[ePos]->size(), "N");


  } else if(group == "FillMetCuts") {
    histAddVal(sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)), "MHT");
    histAddVal(sumptForHt, "HT");  
    histAddVal(sumptForHt + sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)), "Meff");
    histAddVal(theMETVector.Pt(), "Met");
    
  } else if(group == "FillLeadingJet" && goodParts[CUTS::eSusyCom]->size() == 0) {

    if(goodParts[CUTS::eR1stJet]->at(0) != -1) {
      histAddVal(_Jet->smearP.at(goodParts[CUTS::eR1stJet]->at(0)).Pt(), "FirstPt");
      histAddVal(_Jet->smearP.at(goodParts[CUTS::eR1stJet]->at(0)).Eta(), "FirstEta");
    }
    if(goodParts[CUTS::eR2ndJet]->at(0) != -1) {
      histAddVal(_Jet->smearP.at(goodParts[CUTS::eR2ndJet]->at(0)).Pt(), "SecondPt");
      histAddVal(_Jet->smearP.at(goodParts[CUTS::eR2ndJet]->at(0)).Eta(), "SecondEta");
    }

  
  } else if(group == "FillLeadingJet" && goodParts[CUTS::eSusyCom]->size() != 0) {

    TLorentzVector first = _Jet->smearP.at(goodParts[CUTS::eR1stJet]->at(0));
    TLorentzVector second = _Jet->smearP.at(goodParts[CUTS::eR2ndJet]->at(0));
    
    histAddVal(first.Pt(), "FirstPt");
    histAddVal(second.Pt(), "SecondPt");

    histAddVal(first.Eta(), "FirstEta");
    histAddVal(second.Eta(), "SecondEta");
    
    TLorentzVector LeadDiJet = first + second;
    
    histAddVal(LeadDiJet.M(), "Mass"); 
    histAddVal(LeadDiJet.Pt(), "Pt");  
    histAddVal(fabs(first.Eta() - second.Eta()), "DeltaEta"); 
    histAddVal(first.DeltaR(second), "DeltaR");  

    double dphiDijets = absnormPhi(first.Phi() - second.Phi());
    double dphi1 = normPhi(first.Phi() - theMETVector.Phi());
    double dphi2 = normPhi(second.Phi() - theMETVector.Phi());
    double alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histAddVal(dphiDijets, "LeadSublDijetDphi"); 
    histAddVal2(theMETVector.Pt(),dphiDijets, "MetVsDiJetDeltaPhiLeadSubl");
    histAddVal2(fabs(first.Eta()-second.Eta()), dphiDijets, "DeltaEtaVsDeltaPhiLeadSubl");

    histAddVal(absnormPhi(theMETVector.Phi() - LeadDiJet.Phi()), "MetDeltaPhi");



    histAddVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), "R1");
    histAddVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), "R2");
    histAddVal(normPhi(first.Phi() - phiForMht), "Dphi1MHT"); 
    histAddVal(normPhi(second.Phi() - phiForMht), "Dphi2MHT");
    histAddVal(dphi1, "Dphi1");
    histAddVal(dphi2, "Dphi2");
    histAddVal2(dphi1,dphi2, "Dphi1VsDphi2");
    histAddVal(alpha, "Alpha");


    //dijet info
  } else if(group == "FillDiJet") {
    double leaddijetmass = 0;
    double leaddijetpt = 0;
    double leaddijetdeltaR = 0;
    double leaddijetdeltaEta = 0;
    double etaproduct = 0;
    for(vec_iter it=goodParts[CUTS::eDiJet]->begin(); it!=goodParts[CUTS::eDiJet]->end(); it++) {
      int p1 = (*it) / _Jet->smearP.size();
      int p2 = (*it) % _Jet->smearP.size();
      TLorentzVector jet1 = _Jet->smearP.at(p1);
      TLorentzVector jet2 = _Jet->smearP.at(p2);
      TLorentzVector DiJet = jet1 + jet2;
	  
      if(DiJet.M() > leaddijetmass) {
	leaddijetmass = DiJet.M();
	etaproduct = (jet1.Eta() * jet2.Eta() > 0) ? 1 : -1;
      }
      if(DiJet.Pt() > leaddijetpt) leaddijetpt = DiJet.Pt();
      if(fabs(jet1.Eta() - jet2.Eta()) > leaddijetdeltaEta) leaddijetdeltaEta = fabs(jet1.Eta() - jet2.Eta());
      if(jet1.DeltaR(jet2) > leaddijetdeltaR) leaddijetdeltaR = jet1.DeltaR(jet2);

      histAddVal(DiJet.M(), "Mass");
      histAddVal(DiJet.Pt(), "Pt");
      histAddVal(fabs(jet1.Eta() - jet2.Eta()), "DeltaEta");
      histAddVal(absnormPhi(jet1.Phi() - jet2.Phi()), "DeltaPhi");
      histAddVal(jet1.DeltaR(jet2), "DeltaR");
    }


    histAddVal(leaddijetmass, "LargestMass");
    histAddVal(leaddijetpt, "LargestPt");  
    histAddVal(leaddijetdeltaEta, "LargestDeltaEta");
    histAddVal(leaddijetdeltaR, "LargestDeltaR");
    histAddVal(etaproduct, "LargestMassEtaProduct");


    ////diparticle stuff
  } else if(fillInfo[group]->type == FILLER::Dipart) {
    Lepton* lep1 = static_cast<Lepton*>(fillInfo[group]->part);
    Lepton* lep2 = static_cast<Lepton*>(fillInfo[group]->part2);
    CUTS ePos = fillInfo[group]->ePos;
    string digroup = group;
    digroup.erase(0,4);

    TLorentzVector part1;
    TLorentzVector part2;
    
    for(vec_iter it=goodParts[ePos]->begin(); it!=goodParts[ePos]->end(); it++) {
      
      int p1= (*it) / BIG_NUM;
      int p2= (*it) % BIG_NUM;

      part1 = lep1->smearP.at(p1);
      part2 = lep2->smearP.at(p2);
      double add_px = (part1.Pt() * cos(part1.Phi())) + (part2.Pt() * cos(part2.Phi()));
      double add_py = (part1.Pt() * sin(part1.Phi())) + (part2.Pt() * sin(part2.Phi()));
//    TLorentzVector add_p = part1 + part2;
      double mag_addp2 = pow(add_px, 2) + pow(add_py, 2);


      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR"); 
      if(group.find("Di") != string::npos) {
	histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");  
	histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
	histAddVal(TMath::Sqrt(mag_addp2), "Boost");
      } else {
	histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");  
	histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
	histAddVal(TMath::Sqrt(mag_addp2), "Boost");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");
      histAddVal(absnormPhi(part1.Phi() - theMETVector.Phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - theMETVector.Phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - theMETVector.Phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - theMETVector.Phi())), "CosDphi_DeltaPtAndMet");
      double diMass = diParticleMass(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"));
      if(passDiParticleApprox(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"))) {
	histAddVal(diMass, "ReconstructableMass");
      } else {
      	histAddVal(diMass, "NotReconstructableMass");
      }

      double PZeta = getPZeta(part1,part2).first;
      double PZetaVis = getPZeta(part1,part2).second;
      histAddVal(calculateLeptonMetMt(part1), "Part1MetMt");
      histAddVal(calculateLeptonMetMt(part2), "Part2MetMt"); 
      histAddVal(lep2->charge->at(p2) * lep1->charge->at(p1), "OSLS");  
      histAddVal(PZeta, "PZeta"); 
      histAddVal(PZetaVis, "PZetaVis");  
      histAddVal2(PZetaVis,PZeta, "Zeta2D");  
      histAddVal((distats.at(digroup).dmap.at("PZetaCutCoefficient") * PZeta) + (distats.at(digroup).dmap.at("PZetaVisCutCoefficient") * PZetaVis), "Zeta1D");

      if ((goodParts[CUTS::eR1stJet]->at(0) != -1) && (goodParts[CUTS::eR2ndJet]->at(0) != -1)) {
	TLorentzVector TheLeadDiJetVect = _Jet->smearP.at(goodParts[CUTS::eR1stJet]->at(0)) + _Jet->smearP.at(goodParts[CUTS::eR2ndJet]->at(0));

	histAddVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), "Part1DiJetDeltaPhi");
	histAddVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), "Part2DiJetDeltaPhi");
	histAddVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), "DiJetReconstructableMass"); 
      }
      
      if(lep1->type != PType::Tau) {
	histAddVal(isZdecay(part1, *lep1), "Part1IsZdecay"); 
      }
      if(lep2->type != PType::Tau){ 
	histAddVal(isZdecay(part2, *lep2), "Part2IsZdecay"); 
      }
    }
  }
}


void Analyzer::initializePileupInfo(string MCHisto, string DataHisto) {

  TFile *file1 = new TFile((PUSPACE+MCHisto).c_str());
  TH1D* histmc = (TH1D*)file1->FindObjectAny("pileup");
  if(!histmc) throw std::runtime_error("failed to extract histogram");

  TFile* file2 = new TFile((PUSPACE+DataHisto).c_str());
  TH1D* histdata = (TH1D*)file2->FindObjectAny("pileup");
  if(!histdata) throw std::runtime_error("failed to extract histogram");

  double factor = histmc->Integral() / histdata->Integral();
  double value;
  for(int bin=0; bin < 100; bin++) {
    if(histmc->GetBinContent(bin) == 0) value = 1;

    else value = factor*histdata->GetBinContent(bin) / histmc->GetBinContent(bin);
    hPU[bin] = value;
  }

  file1->Close();
  file2->Close();
  
}

double normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}

double absnormPhi(double phi) {
  return abs(normPhi(phi));
}

bool Analyzer::partPassBoth(string name) {
  int diffsize = 0;
  CUTS ePart1;
  CUTS ePart2;
  if(name == "Muon1Muon2") {
    ePart1 = CUTS::eRMuon1;
    ePart2 = CUTS::eRMuon2;
    diffsize = _Muon->smearP.size();      
  } else if(name == "Electron1Electron2") {
    ePart1 = CUTS::eRElec1;
    ePart2 = CUTS::eRElec2;
    diffsize = _Electron->smearP.size();
  } else if(name == "Tau1Tau2") {
    ePart1 = CUTS::eRTau1;
    ePart2 = CUTS::eRTau2;
    diffsize = _Tau->smearP.size();
  }
  vector<int>* part1 = goodParts[ePart1];
  vector<int>* part2 = goodParts[ePart2];
  vector<int> diff(diffsize);

  vector<int>::iterator it = set_symmetric_difference(part1->begin(), part1->end(), part2->begin(), part2->end(), diff.begin());
  diff.resize(it - diff.begin());

  return (diff.size() == 0);
}


bool CRTester::test(Analyzer* analyzer) {

  bool pass = false;

  if(info->type == FILLER::Single) {
    if(variable == "PassBoth") return analyzer->partPassBoth(partName);

    for(auto index: *analyzer->getList(info->ePos)) {
  
      TLorentzVector part = info->part->smearP.at(index);
      if(variable == "Eta") pass = pass && part.Eta() > cutVal;
      else if(variable == "Pt") pass = pass && part.Pt() > cutVal;
      else if(variable == "Energy") pass = pass && part.Energy() > cutVal;

      //	else if(variable == "Zdecay") pass = pass && analyzer.isZdecay(part, static_cast<Lepton>(*info->part)) > cutVal;
      //else if(variable == "MetPhi") pass = pass && part.Eta() > cutval;
    }	  
  } else if(info->type == FILLER::Dipart) {
    for(auto index: *analyzer->getList(info->ePos)) {

      TLorentzVector part1 = info->part->smearP.at(index / BIG_NUM);
      TLorentzVector part2 = info->part2->smearP.at(index % BIG_NUM);
      if(variable == "DeltaR") pass = pass && part1.DeltaR(part2) > cutVal;
      else if(variable == "DeltaPtDivSumPt") pass = pass && ((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()) > cutVal);
      else if(variable == "DeltaPt") pass = pass && ((part1.Pt() - part2.Pt()) > cutVal);
      else if(variable == "Zeta") pass = pass && (analyzer->getZeta(part1, part2, partName) > cutVal);
      //else if(variable == "OSLS") pass = pass && 
      else if(variable == "CosDphi") pass = pass && absnormPhi( part1.Phi() - part2.Phi()) > cutVal;
      else if(variable == "Mass") pass = pass && analyzer->getMass(part1, part2, partName) > cutVal;
      else if(variable == "DeltaEta") pass = pass && (abs(part1.Eta() - part2.Eta()) > cutVal);
      else if(variable == "DeltaPhi") pass = pass && (abs(part1.Phi() - part2.Phi()) > cutVal);
      else if(variable == "OSEta") pass = pass && (part1.Eta() * part2.Eta() > cutVal);
    }   
  } else {

    if(variable == "Met") pass = (analyzer->getMet() > cutVal);
    else if(variable == "Ht") pass = (analyzer->getHT() > cutVal);
    else if(variable == "Mht") pass = (analyzer->getMHT() > cutVal);
  }


  return pass;
}
