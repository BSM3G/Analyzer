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
#define histAddVal2(val1, val2, name) ihisto.addVal(val1, val2, group, max, name, wgt)
#define histAddVal(val, name) ihisto.addVal(val, group, max, name, wgt)
#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

typedef vector<int>::iterator vec_iter;



//////////////////////////////////////////////////////////////////
///////////////////CONSTANTS DEFINITONS///////////////////////////
//////////////////////////////////////////////////////////////////

//Filespace that has all of the .in files
const string PUSPACE = "Pileup/";


//////////PUBLIC FUNCTIONS////////////////////

const vector<CUTS> Analyzer::genCuts = {
  CUTS::eGTau, CUTS::eNuTau, CUTS::eGTop,
  CUTS::eGElec, CUTS::eGMuon, CUTS::eGZ,
  CUTS::eGW, CUTS::eGHiggs, CUTS::eGJet
};

const vector<CUTS> Analyzer::jetCuts = {
  CUTS::eRJet1,  CUTS::eRJet2,   CUTS::eRCenJet,
  CUTS::eR1stJet, CUTS::eR2ndJet, CUTS::eRBJet
};

const vector<CUTS> Analyzer::nonParticleCuts = {
  CUTS::eRVertex,CUTS::eRTrig1, CUTS::eRTrig2,
};

const unordered_map<string, CUTS> Analyzer::cut_num = {
  {"NGenTau", CUTS::eGTau},                             {"NGenTop", CUTS::eGTop},
  {"NGenElectron", CUTS::eGElec},                       {"NGenMuon", CUTS::eGMuon},
  {"NGenZ", CUTS::eGZ},                                 {"NGenW", CUTS::eGW},
  {"NGenHiggs", CUTS::eGHiggs},                         {"NGenJet", CUTS::eGJet},
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
  {"NLeadJetCombinations", CUTS::eSusyCom},             {"METCut", CUTS::eMET},
  {"NRecoWJet", CUTS::eRWjet},                          {"NRecoVertex", CUTS::eRVertex}
};


//////////////////////////////////////////////////////
//////////////////PUBLIC FUNCTIONS////////////////////
//////////////////////////////////////////////////////

///Constructor
Analyzer::Analyzer(vector<string> infiles, string outfile, bool setCR, string configFolder) : goodParts(getArray()), genName_regex(".*([A-Z][^[:space:]]+)") {
  cout << "setup start" << endl;

  BOOM= new TChain("TNT/BOOM");
  infoFile=0;

  for( string infile: infiles){
    BOOM->AddFile(infile.c_str());
  }


  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "TOTAL EVENTS: " << nentries << std::endl;

  srand(0);

  for(int i=0; i < nTrigReq; i++) {
    vector<int>* tmpi = new vector<int>();
    vector<string>* tmps = new vector<string>();
    trigPlace[i] = tmpi;
    trigName[i] = tmps;
  }

  filespace=configFolder;
  filespace+="/";

  setupGeneral();

  reader.load(calib, BTagEntry::FLAV_B, "comb");

  isData = distats["Run"].bfind("isData");

  CalculatePUSystematics = distats["Run"].bfind("CalculatePUSystematics");
  initializePileupInfo(distats["Run"].smap.at("MCHistos"), distats["Run"].smap.at("DataHistos"),distats["Run"].smap.at("DataPUHistName"),distats["Run"].smap.at("MCPUHistName"));
  syst_names.push_back("orig");
  unordered_map<CUTS, vector<int>*, EnumHash> tmp;
  syst_parts.push_back(tmp);
  if(!isData && distats["Systematics"].bfind("useSystematics")) {
    for(auto systname : distats["Systematics"].bset) {
      if( systname == "useSystematics")
        doSystematics= true;
      else {
        syst_names.push_back(systname);
        syst_parts.push_back(getArray());
      }
    }
  }else {
    doSystematics=false;
  }
  cout << doSystematics << endl;
  _Electron = new Electron(BOOM, filespace + "Electron_info.in", syst_names);
  _Muon = new Muon(BOOM, filespace + "Muon_info.in", syst_names);
  _Tau = new Taus(BOOM, filespace + "Tau_info.in", syst_names);
  _Jet = new Jet(BOOM, filespace + "Jet_info.in", syst_names);
  _FatJet = new FatJet(BOOM, filespace + "FatJet_info.in", syst_names);
  _MET  = new Met(BOOM, "Met_type1PF" , syst_names);

  if(!isData) {
    _Gen = new Generated(BOOM, filespace + "Gen_info.in", syst_names);
    allParticles= {_Gen,_Electron,_Muon,_Tau,_Jet,_FatJet};
  } else {
    allParticles= {_Electron,_Muon,_Tau,_Jet,_FatJet};
  }


  vector<string> cr_variables;
  if(setCR) {
    char buf[64];
    read_info(filespace + "Control_Regions.in");
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
      if(distats["Control_Region"].smap.find("Unblind") != distats["Control_Region"].smap.end()) {

        blinded = distats["Control_Region"].smap["Unblind"] == "false";
        cout << "we have " << blinded << endl;
      }
    }
  }
  //we update the root file if it exist so now we have to delete it:
  std::remove(outfile.c_str()); // delete file
  histo = Histogramer(1, filespace+"Hist_entries.in", filespace+"Cuts.in", outfile, isData, cr_variables);
  if(doSystematics)
    //cout << "syst_names: " << syst_names.size() << endl;
    syst_histo=Histogramer(1, filespace+"Hist_syst_entries.in", filespace+"Cuts.in", outfile, isData, cr_variables,syst_names);
  systematics = Systematics(distats);
  jetScaleRes = JetScaleResolution("Pileup/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt", "",  "Pileup/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt", "Pileup/Spring16_25nsV6_MC_SF_AK4PFchs.txt");





  ///this can be done nicer
  //put the variables that you use here:
  //zBoostTree["treeweight"] = 0;
  //zBoostTree["mu1pt"] = 0;
  //zBoostTree["mu2pt"] = 0;
  //zBoostTree["addpxZ"] = 0;
  //zBoostTree["addpyZ"] = 0;
  //zBoostTree["utpx"] = 0;
  //zBoostTree["utpy"] = 0;
  //zBoostTree["utmagn"] = 0;
  //zBoostTree["utpar"] = 0;
  //zBoostTree["utper"] = 0;
  //zBoostTree["boostedz"] = 0;
  //zBoostTree["FLJetPt"] = 0;

  //histo.createTree(&zBoostTree,"zboost");


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
  // check if we need to make gen level cuts to cross clean the samples:

  isVSample = false;
  if(infiles[0].find("DY") != string::npos){
    isVSample = true;
    if(infiles[0].find("DYJetsToLL_M-50_HT-") != string::npos){
      //gen_selection["DY_noMass_gt_100"]=true;
      gen_selection["DY_noMass_gt_200"]=true;
    //get the DY1Jet DY2Jet ...
    }else if(infiles[0].find("JetsToLL_TuneCUETP8M1_13TeV") != string::npos){
      gen_selection["DY_noMass_gt_200"]=true;
    }else{
      //set it to false!!
      gen_selection["DY_noMass_gt_200"]=false;
    }
    if(infiles[0].find("DYJetsToLL_M-50_TuneCUETP8M1_13TeV") != string::npos){
      gen_selection["DY_noMass_gt_200"]=true;
    }else{
      //set it to false!!
      gen_selection["DY_noMass_gt_100"]=false;
    }
  }else{
    //set it to false!!
    gen_selection["DY_noMass_gt_200"]=false;
    gen_selection["DY_noMass_gt_100"]=false;
  }

  if(infiles[0].find("WJets") != string::npos){
    isVSample = true;
  }

  for(auto iselect : gen_selection){
    if(iselect.second){
      cout<<"Waning: The selection "<< iselect.first<< " is active!"<<endl;
    }
  }

  setCutNeeds();

  std::cout << "setup complete" << std::endl << endl;
  start = std::chrono::system_clock::now();
}

unordered_map<CUTS, vector<int>*, EnumHash> Analyzer::getArray() {
  unordered_map<CUTS, vector<int>*, EnumHash> rmap;
  for(auto e: Enum<CUTS>()) {
    rmap[e] = new vector<int>();
  }
  return rmap;
}



void Analyzer::create_fillInfo() {

  fillInfo["FillLeadingJet"] = new FillVals(CUTS::eSusyCom, FILLER::Dipart, _Jet, _Jet);
  fillInfo["FillGen"] =        new FillVals(CUTS::eGen, FILLER::Single, _Gen);
  fillInfo["FillTau1"] =       new FillVals(CUTS::eRTau1, FILLER::Single, _Tau);
  fillInfo["FillTau2"] =       new FillVals(CUTS::eRTau2, FILLER::Single, _Tau);
  fillInfo["FillMuon1"] =      new FillVals(CUTS::eRMuon1, FILLER::Single, _Muon);
  fillInfo["FillMuon2"] =      new FillVals(CUTS::eRMuon2, FILLER::Single, _Muon);
  fillInfo["FillElectron1"] =  new FillVals(CUTS::eRElec1, FILLER::Single, _Electron);
  fillInfo["FillElectron2"] =  new FillVals(CUTS::eRElec2, FILLER::Single, _Electron);

  fillInfo["FillJet1"] =       new FillVals(CUTS::eRJet1, FILLER::Single, _Jet);
  fillInfo["FillJet2"] =       new FillVals(CUTS::eRJet2, FILLER::Single, _Jet);
  fillInfo["FillBJet"] =       new FillVals(CUTS::eRBJet, FILLER::Single, _Jet);
  fillInfo["FillCentralJet"] = new FillVals(CUTS::eRCenJet, FILLER::Single, _Jet);
  fillInfo["FillWJet"] =     new FillVals(CUTS::eRWjet, FILLER::Single, _FatJet);

  fillInfo["FillDiElectron"] = new FillVals(CUTS::eDiElec, FILLER::Dipart, _Electron, _Electron);
  fillInfo["FillDiMuon"] =     new FillVals(CUTS::eDiMuon, FILLER::Dipart, _Muon, _Muon);
  fillInfo["FillDiTau"] =      new FillVals(CUTS::eDiTau, FILLER::Dipart, _Tau, _Tau);
  fillInfo["FillMetCuts"] =    new FillVals();
  fillInfo["FillDiJet"] =      new FillVals(CUTS::eDiJet, FILLER::Dipart, _Jet, _Jet);

  fillInfo["FillMuon1Tau1"] =       new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon1Tau2"] =       new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau1"] =       new FillVals(CUTS::eMuon2Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau2"] =       new FillVals(CUTS::eMuon2Tau2, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillElectron1Tau1"] =   new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron1Tau2"] =   new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau1"] =   new FillVals(CUTS::eElec2Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau2"] =   new FillVals(CUTS::eElec2Tau2, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillMuon1Electron1"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon1Electron2"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron1"] =  new FillVals(CUTS::eMuon2Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron2"] =  new FillVals(CUTS::eMuon2Elec2, FILLER::Dipart, _Muon, _Electron);

  //////I hate this solution so much.  Its terrible
  fillInfo["FillElectron1Electron2"] =     new FillVals(CUTS::eDiElec, FILLER::Single, _Electron, _Electron);
  fillInfo["FillMuon1Muon2"] =             new FillVals(CUTS::eDiMuon, FILLER::Single, _Muon, _Muon);
  fillInfo["FillTau1Tau2"] =               new FillVals(CUTS::eDiTau, FILLER::Single, _Tau, _Tau);



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
  clear_values();
  delete BOOM;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
  if(!isData) delete _Gen;

  for(auto pair: fillInfo) {
    delete pair.second;
    pair.second=nullptr;
  }

  for(auto e: Enum<CUTS>()) {
    delete goodParts[e];
    goodParts[e]=nullptr;
  }
  //for(auto &it: syst_parts) {
    //for(auto e: Enum<CUTS>()) {
      //if( it[e] != nullptr) {
      //if(it.find(e) != it.end()){
        //delete it[e];
        //it[e]=nullptr;
      //}
      //}
    //}
  //}
  for(auto it: testVec){
    delete it;
    it=nullptr;
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
  //faster!!
  for(auto &it: syst_parts) {
    if (it.size() == 0) continue;
    for(auto e: Enum<CUTS>()) {
      it[e]->clear();
    }
  }
  if(infoFile!=BOOM->GetFile()){
    cout<<"New file!"<<endl;
    infoFile=BOOM->GetFile();
  }
  if(version==1 && infoFile!=BOOM->GetFile()){
    cout<<"New file! Will get the trigger info."<<endl;
    infoFile=BOOM->GetFile();
    BAAM= (TTree*) infoFile->Get("TNT/BAAM");
    initializeTrigger();
    infoFile->Close();
  }


  leadIndex=-1;
  maxCut = 0;
}

///Function that does most of the work.  Calculates the number of each particle
void Analyzer::preprocess(int event) {
  BOOM->GetEntry(event);
  for(Particle* ipart: allParticles){
    ipart->init();
  }
  _MET->init();

  active_part = &goodParts;
  if(!select_mc_background()){
    //we will put nothing in good particles
    clear_values();
    return;
  }

  pu_weight = (!isData && CalculatePUSystematics) ? hPU[(int)(nTruePU+1)] : 1.0;
  
  // SET NUMBER OF GEN PARTICLES
  if(!isData){
    _Gen->setOrigReco();
    getGoodGen(_Gen->pstats["Gen"]);
    getGoodTauNu();
    //---MY APPLICATION OF BOOST WEIGHTS-------------------------------------------
  //   bool checkw = checkforW();  //new7.27.17 See if there's a W with the function checkforW.
  //   bool checkz = checkforZ(); //new9.15.17 See if there's a Z with the function checkforZ.
  //   if(checkw){
  //     double boostw = getWBoostWeight(); //Grab the p_T of the W.
  //     if(boostw > 0 && boostw <= 50) {pu_weight *= 1.1192;}  //new7.27.17
  //     else if (boostw > 50 && boostw <= 100) {pu_weight *= 1.1034;}   //new7.27.17
  //     else if (boostw > 100 && boostw <= 150) {pu_weight *= 1.0675;}  //new7.27.17
  //     else if (boostw > 150 && boostw <= 200) {pu_weight *= 1.0637;}  //new7.27.17
  //     else if (boostw > 200 && boostw <= 300) {pu_weight *= 1.0242;}  //new7.27.17
  //     else if (boostw > 300 && boostw <= 400) {pu_weight *= 0.9453;}  //new7.27.17
  //     else if (boostw > 400 && boostw <= 600) {pu_weight *= 0.8579;}  //new7.27.17
  //     else if (boostw >= 600) {pu_weight *= 0.7822;}  //new7.27.17
  //     else {pu_weight *= 1;}  //new7.27.17
  //   }
  //   else if (checkz){
  //     double boostz = getZBoostWeight(); //Grab the p_T of the Z.
  //     if(boostz > 0 && boostz <= 50) {pu_weight *= 1.1192;}  //new9.15.17
  //     else if (boostz > 50 && boostz <= 100) {pu_weight *= 1.1034;}   //new9.15.17                                                                                 
  //     else if (boostz > 100 && boostz <= 150) {pu_weight *= 1.0675;}  //new9.15.17                                                                                 
  //     else if (boostz > 150 && boostz <= 200) {pu_weight *= 1.0637;}  //new9.15.17                                                                                 
  //     else if (boostz > 200 && boostz <= 300) {pu_weight *= 1.0242;}  //new9.15.17                                                                                 
  //     else if (boostz > 300 && boostz <= 400) {pu_weight *= 0.9453;}  //new9.15.17                                                                                 
  //     else if (boostz > 400 && boostz <= 600) {pu_weight *= 0.8579;}  //new9.15.17                                                                                 
  //     else if (boostz >= 600) {pu_weight *= 0.7822;}  //new9.15.17           
  //     else {pu_weight *= 1;}  //new9.15.17  
  //   }
  //   else {pu_weight *= 1;}  //new7.27.17
  }
    //---END MY APPLICATION OF BOOST WEIGHTS---------------------------------------
  //////Triggers and Vertices
  active_part->at(CUTS::eRVertex)->resize(bestVertices);
  TriggerCuts(*(trigPlace[0]), *(trigName[0]), CUTS::eRTrig1);
  TriggerCuts(*(trigPlace[1]), *(trigName[1]), CUTS::eRTrig2);

  ////check update met is ok  
  for(size_t i=0; i < syst_names.size(); i++) {
     //////Smearing
    smearLepton(*_Electron, CUTS::eGElec, _Electron->pstats["Smear"], distats["Electron_systematics"], i);
    smearLepton(*_Muon, CUTS::eGMuon, _Muon->pstats["Smear"], distats["Muon_systematics"], i);
    smearLepton(*_Tau, CUTS::eGTau, _Tau->pstats["Smear"], distats["Tau_systematics"], i);

    smearJet(*_Jet,CUTS::eGJet,_Jet->pstats["Smear"], i);
    smearJet(*_FatJet,CUTS::eGJet,_FatJet->pstats["Smear"], i);
    updateMet(i);
    
  }
  
  for(size_t i=0; i < syst_names.size(); i++) {
    string systname = syst_names.at(i);
    for( auto part: allParticles) part->setCurrentP(i);
    _MET->setCurrentP(i);
    getGoodParticles(i);
  }
  active_part = &goodParts;
  if( event < 10 || ( event < 100 && event % 10 == 0 ) ||
    ( event < 1000 && event % 100 == 0 ) ||
    ( event < 10000 && event % 1000 == 0 ) ||
    ( event >= 10000 && event % 10000 == 0 ) ) {
       cout << setprecision(2)<<event << " Events analyzed "<< static_cast<double>(event)/nentries*100. <<"% done"<<endl;
       cout << fixed;
  }
}

void Analyzer::getGoodParticles(int syst){

  string systname=syst_names.at(syst);
  if(syst == 0) active_part = &goodParts;
  else active_part=&syst_parts.at(syst);
    //    syst=syst_names[syst];


  // // SET NUMBER OF RECO PARTICLES
  // // MUST BE IN ORDER: Muon/Electron, Tau, Jet
  getGoodRecoLeptons(*_Electron, CUTS::eRElec1, CUTS::eGElec, _Electron->pstats["Elec1"],syst);
  getGoodRecoLeptons(*_Electron, CUTS::eRElec2, CUTS::eGElec, _Electron->pstats["Elec2"],syst);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon1, CUTS::eGMuon, _Muon->pstats["Muon1"],syst);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon2, CUTS::eGMuon, _Muon->pstats["Muon2"],syst);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau1, CUTS::eGTau, _Tau->pstats["Tau1"],syst);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau2, CUTS::eGTau, _Tau->pstats["Tau2"],syst);

  getGoodRecoJets(CUTS::eRBJet, _Jet->pstats["BJet"],syst);
  getGoodRecoJets(CUTS::eRJet1, _Jet->pstats["Jet1"],syst);
  getGoodRecoJets(CUTS::eRJet2, _Jet->pstats["Jet2"],syst);
  getGoodRecoJets(CUTS::eRCenJet, _Jet->pstats["CentralJet"],syst);
  getGoodRecoJets(CUTS::eR1stJet, _Jet->pstats["FirstLeadingJet"],syst);
  getGoodRecoJets(CUTS::eR2ndJet, _Jet->pstats["SecondLeadingJet"],syst);

  getGoodRecoFatJets(CUTS::eRWjet, _FatJet->pstats["Wjet"],syst);
  //  treatMuons_Met(systname);


  ///VBF Susy cut on leadin jets
  VBFTopologyCut(distats["VBFSUSY"],syst);

  /////lepton lepton topology cuts
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1,CUTS::eRTau1, CUTS::eElec1Tau1, distats["Electron1Tau1"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau1, CUTS::eElec2Tau1, distats["Electron2Tau1"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1, CUTS::eRTau2, CUTS::eElec1Tau2, distats["Electron1Tau2"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau2, CUTS::eElec2Tau2, distats["Electron2Tau2"],syst);

  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau1, CUTS::eMuon1Tau1, distats["Muon1Tau1"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau2, CUTS::eMuon1Tau2, distats["Muon1Tau2"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau1, CUTS::eMuon2Tau1, distats["Muon2Tau1"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau2, CUTS::eMuon2Tau2, distats["Muon2Tau2"],syst);

  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec1, CUTS::eMuon1Elec1, distats["Muon1Electron1"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec2, CUTS::eMuon1Elec2, distats["Muon1Electron2"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec1, CUTS::eMuon2Elec1, distats["Muon2Electron1"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec2, CUTS::eMuon2Elec2, distats["Muon2Electron2"],syst);

  ////DIlepton topology cuts
  getGoodLeptonCombos(*_Tau, *_Tau, CUTS::eRTau1, CUTS::eRTau2, CUTS::eDiTau, distats["DiTau"],syst);
  getGoodLeptonCombos(*_Electron, *_Electron, CUTS::eRElec1, CUTS::eRElec2, CUTS::eDiElec, distats["DiElectron"],syst);
  getGoodLeptonCombos(*_Muon, *_Muon, CUTS::eRMuon1, CUTS::eRMuon2, CUTS::eDiMuon, distats["DiMuon"],syst);

  ////Dijet cuts
  getGoodDiJets(distats["DiJet"],syst);

}

////Reads cuts from Cuts.in file and see if the event has enough particles
bool Analyzer::fillCuts(bool fillCounter) {
  const unordered_map<string,pair<int,int> >* cut_info = histo.get_cuts();
  const vector<string>* cut_order = histo.get_cutorder();
 
  bool prevTrue = true;
  
  maxCut=0;
  //  cout << active_part << endl;;

  for(size_t i = 0; i < cut_order->size(); i++) {
    string cut = cut_order->at(i);
    if(isData && cut.find("Gen") != string::npos) continue;
    int min= cut_info->at(cut).first;
    int max= cut_info->at(cut).second;
    int nparticles = active_part->at(cut_num.at(cut))->size();
    //if(!fillCounter) cout << cut << ": " << nparticles << " (" << min << ", " << max << ")" <<endl;
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      if((cut_num.at(cut) == CUTS::eR1stJet || cut_num.at(cut) == CUTS::eR2ndJet) && active_part->at(cut_num.at(cut))->at(0) == -1 ) {
        //cout<<"here   "<<endl;
        prevTrue = false;
        continue;  ////dirty dirty hack
      }
      if(fillCounter && crbins == 1) {
        cuts_per[i]++;
        cuts_cumul[i] += (prevTrue) ? 1 : 0;
        maxCut += (prevTrue) ? 1 : 0;
      }
    }else {
      //cout<<"here 2  "<<endl;
      prevTrue = false;
    }
  }
  
  if(crbins != 1) {
    if(!prevTrue) {
      maxCut = -1;
      return prevTrue;
    }

    int factor = crbins;
    for(auto tester: testVec) {
      factor /= 2;
      /////get variable value from maper.first.
      if(tester->test(this)) { ///pass cut
        maxCut += factor;
      }
    }
    if(isData && blinded && maxCut == SignalRegion) return false;
    cuts_per[maxCut]++;
  }

  
  return prevTrue;
}



///Prints the number of events that passed each cut per event and cumulatively
//done at the end of the analysis
void Analyzer::printCuts() {
  vector<string> cut_order;
  if(crbins > 1) cut_order = *(histo.get_folders());
  else cut_order = *(histo.get_cutorder());
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  double run_time_real=elapsed_seconds.count();


  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(3);
  cout << "\n";
  cout << "Selection Efficiency " << "\n";
  cout << "Total events: " << nentries << "\n";
  cout << "\n";
  cout << "Run Time (real): " <<run_time_real <<" s\n";
  cout << "Time per 1k Events (real): " << run_time_real/(nentries/1000) <<" s\n";
  cout << "Events/s: " << static_cast<double>(nentries)/(run_time_real) <<" 1/s (real) \n";
  cout << "                        Name                  Indiv.";
  if(crbins == 1) cout << "            Cumulative";
  cout << endl << "---------------------------------------------------------------------------\n";
  for(size_t i = 0; i < cut_order.size(); i++) {
    cout << setw(28) << cut_order.at(i) << "    ";
    if(isData && cut_order.at(i).find("Gen") != string::npos) cout << "Skipped" << endl;
    else if(crbins != 1 && blinded && i == (size_t)SignalRegion) cout << "Blinded Signal Region" << endl;
    else {
      cout << setw(10) << cuts_per.at(i) << "  ( " << setw(5) << ((float)cuts_per.at(i)) / nentries << ") ";
      if(crbins == 1) cout << setw(12) << cuts_cumul.at(i) << "  ( " << setw(5) << ((float)cuts_cumul.at(i)) / nentries << ") ";

      cout << endl;
    }
  }
  cout << "---------------------------------------------------------------------------\n";

  //write all the histograms
  //attention this is not the fill_histogram method from the Analyser
  cout << "--------\n" << endl;
  histo.fill_histogram();
  if(doSystematics)
    syst_histo.fill_histogram();

}

/////////////PRIVATE FUNCTIONS////////////////


bool Analyzer::select_mc_background(){
  //will return true if Z* mass is smaller than 200GeV
  if(gen_selection["DY_noMass_gt_200"]){
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      if(abs(_Gen->pdg_id->at(i))==11 or abs(_Gen->pdg_id->at(i))==13 or abs(_Gen->pdg_id->at(i))==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          return (lep1+lep2).M()<200;
        }else{
          lep1= _Gen->p4(i);
        }
      }
    }
  }
  //will return true if Z* mass is smaller than 200GeV
  if(gen_selection["DY_noMass_gt_100"]){
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      if(abs(_Gen->pdg_id->at(i))==11 or abs(_Gen->pdg_id->at(i))==13 or abs(_Gen->pdg_id->at(i))==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          return (lep1+lep2).M()<100;
        }else{
          lep1= _Gen->p4(i);
        }
      }
    }
  }
  //cout<<"Something is rotten in the state of Denmark."<<endl;
  //cout<<"could not find gen selection particle"<<endl;
  return true;
}


double Analyzer::getTauDataMCScaleFactor(int updown){
  double sf=1.;
  //for(size_t i=0; i<_Tau->size();i++){
  for(auto i : *active_part->at(CUTS::eRTau1)){
    if(matchTauToGen(_Tau->p4(i),0.4)!=TLorentzVector()){
      
      if(updown==-1) sf*=  _Tau->pstats["Smear"].dmap.at("TauSF") * (1.-(0.35*_Tau->pt(i)/1000.0));
      else if(updown==0) sf*=  _Tau->pstats["Smear"].dmap.at("TauSF");
      else if(updown==1) sf*=  _Tau->pstats["Smear"].dmap.at("TauSF") * (1.+(0.05*_Tau->pt(i)/1000.0));
    }
  }
  return sf;
}

///Calculates met from values from each file plus smearing and treating muons as neutrinos
void Analyzer::updateMet(int syst) {
  _MET->update(distats["Run"], *_Jet,  syst);

  /////MET CUTS

  if(!passCutRange(_MET->pt(), distats["Run"].pmap.at("MetCut"))) return;
  if(distats["Run"].bfind("DiscrByHT") && _MET->HT() < distats["Run"].dmap.at("HtCut")) return;

  active_part->at(CUTS::eMET)->push_back(1);
}

///////////////////////////////////////////////
////////removed for teh time being/////////////
///////////////////////////////////////////////

// void Analyzer::treatMuons_Met(string syst) {

//   //syst not implemented for muon as tau or neutrino yet
//   if( syst!="orig" or !( distats["Run"].bfind("TreatMuonsAsNeutrinos") || distats["Run"].bfind("TreatMuonsAsTaus")) ){
//     return;
//   }

//   //  Neutrino update before calculation
//   _MET->addP4Syst(_MET->p4(),"muMET");
//   _MET->systdeltaMEx["muMET"]=0;
//   _MET->systdeltaMEy["muMET"]=0;

//   if(distats["Run"].bfind("TreatMuonsAsNeutrinos")) {
//     for(auto it : *active_part->at(CUTS::eRMuon1)) {
//       if(find(active_part->at(CUTS::eRMuon2)->begin(), active_part->at(CUTS::eRMuon2)->end(), it) != active_part->at(CUTS::eRMuon2)->end() ) continue;
//       _MET->systdeltaMEx["muMET"] += _Muon->p4(it).Px();
//       _MET->systdeltaMEy["muMET"] += _Muon->p4(it).Py();
//     }
//     for(auto it : *active_part->at(CUTS::eRMuon2)) {
//       _MET->systdeltaMEx["muMET"] += _Muon->p4(it).Px();
//       _MET->systdeltaMEy["muMET"] += _Muon->p4(it).Py();
//     }
//   }
//   // else if(distats["Run"].bmap.at("TreatMuonsAsTaus")) {

//   //   if(active_part->at(CUTS::eRMuon1)->size() == 1) {

//   //     int muon = (int)active_part->at(CUTS::eRMuon1)->at(0);

//   //     double rand1 = 1;//Tau_HFrac->GetRandom();
//   //     double rand2 = 0;//Tau_Resol->GetRandom();

//   //     double ETau_Pt = _Muon->p4(muon).Pt()*rand1*(rand2+1.0);
//   //     double ETau_Eta = _Muon->p4(muon).Eta();
//   //     double ETau_Phi=normPhi(_Muon->p4(muon).Phi());//+DeltaNu_Phi->GetRandom());
//   //     double ETau_Energy = 0.;


//   //     // double theta = 2.0*TMath::ATan2(1.0,TMath::Exp(_Muon->p4(muon).Eta()));
//   //     // double sin_theta = TMath::Sin(theta);
//   //     // double P_tau = ETau_Pt/sin_theta;

//   //     // //ETau_Energy = sqrt(pow(P_tau, 2) + pow(1.77699, 2));
//   //     // ETau_Energy = sqrt( pow(1.77699, 2) + pow(ETau_Pt, 2) + pow(_Muon->p4(muon).Pz(), 2));

//   //     /*if(ETau_Pt <= 15.0){
//   //       while(ETau_Pt<=15.0){
//   //       rand1 = Tau_HFrac->GetRandom();
//   //       rand2 = Tau_Resol->GetRandom();
//   //       ETau_Pt = _Muon->p4(muon).Pt()*rand1*(rand2+1.0);
//   //       ENu_Pt = _Muon->p4(muon).Pt()-ETau_Pt;
//   //       }
//   //     }
//   //     */

//   //     TLorentzVector Emu_Tau;
//   //     Emu_Tau.SetPtEtaPhiE(ETau_Pt, ETau_Eta, ETau_Phi, ETau_Energy);
//   //     _Muon->cur_P->clear();

//   //     if (ETau_Pt >= _Muon->pstats["Muon1"].pmap.at("PtCut").first ){
//   //       _Muon->cur_P->push_back(Emu_Tau);
//   //       _MET->systdeltaMEy["muMET"] += (_Muon->p4(muon).Px()-Emu_Tau.Px());
//   //       _MET->systdeltaMEy["muMET"] += (_Muon->p4(muon).Py()-Emu_Tau.Py());

//   //     }
//   //   }
//   //}
//   // recalculate MET
//   //  _MET->update("muMET");

//   /////MET CUTS
//   active_part->at(CUTS::eMET)->clear();

//   if(passCutRange(_MET->pt(), distats["Run"].pmap.at("MetCut"))) {
//     active_part->at(CUTS::eMET)->push_back(1);
//   }
// }


/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral() {
  
  genMaper = {
    {5, new GenFill(2, CUTS::eGJet)},     {6,  new GenFill(2, CUTS::eGTop)},
    {11, new GenFill(1, CUTS::eGElec)},   {13, new GenFill(1, CUTS::eGMuon)},
    {15, new GenFill(2, CUTS::eGTau)},    {23, new GenFill(62, CUTS::eGZ)},
    {24, new GenFill(62, CUTS::eGW)},      {25, new GenFill(2, CUTS::eGHiggs)}
};

  SetBranch("nTruePUInteractions", nTruePU);
  SetBranch("bestVertices", bestVertices);
  SetBranch("weightevt", gen_weight);
  //SetBranch("rho", rho);

  read_info(filespace + "ElectronTau_info.in");
  read_info(filespace + "MuonTau_info.in");
  read_info(filespace + "MuonElectron_info.in");
  read_info(filespace + "DiParticle_info.in");
  read_info(filespace + "VBFCuts_info.in");
  read_info(filespace + "Run_info.in");
  read_info(filespace + "Systematics_info.in");

  if( BOOM->GetListOfBranches()->FindObject("Trigger_names") ==0){
    SetBranch("Trigger_decision", Trigger_decisionV1);
    infoFile=BOOM->GetFile();
    BAAM= (TTree*) infoFile->Get("TNT/BAAM");
    initializeTrigger();
    infoFile->Close();
    version=1;
  }else{
    SetBranch("Trigger_names", Trigger_names);
    SetBranch("Trigger_decision", Trigger_decision);
    BOOM->GetEntry(0);
    for(int i = 0; i < nTrigReq; i++) {
      for(int j = 0; j < (int)trigName[i]->size(); j++) {
        for(int k = 0; k < (int)Trigger_names->size(); k++) {
          if(Trigger_names->at(k).find(trigName[i]->at(j)) != string::npos) {
            // structure: i tigger 1 or 2 | j  name of trigger in trigger one or two
            trigPlace[i]->at(j) = k;
            break;
          }
        }
      }
    }
    BOOM->SetBranchStatus("Trigger_names", 0);
  }
}


//get the correct trigger position:
void Analyzer::initializeTrigger() {
  BAAM->SetBranchStatus("triggernames", 1);
  BAAM->SetBranchAddress("triggernames", &Trigger_names);

  BAAM->GetEntry(0);
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
  BAAM->SetBranchStatus("triggernames", 0);
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
      if(group.compare("Control_Region") !=0 ){
        if(stemp[1] == "1" || stemp[1] == "true") distats[group].bset.push_back(stemp[0]);
        else if(*p) distats[group].smap[stemp[0]] = stemp[1];
        else  distats[group].dmap[stemp[0]]=stod(stemp[1]);
      }else{
        if(*p) distats[group].smap[stemp[0]] = stemp[1];
        else  distats[group].dmap[stemp[0]]=stod(stemp[1]);
      }

    } else if(stemp.size() == 3) distats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
  }
  info_file.close();
}


// This code works pretty much (at least in my tests), but dagnabit, its ugly.  They all can't be winners, at least now...
void Analyzer::setCutNeeds() {


  for(auto it: *histo.get_groups()) {
    if(fillInfo[it]->type == FILLER::None) continue;
    neededCuts.loadCuts(fillInfo[it]->ePos);
  }
  for(auto it : *histo.get_cutorder()) {
    try{
      neededCuts.loadCuts(cut_num.at(it));
    }catch(...){
      cout<<"The following cut is strange: "<<it<<endl;
      exit(2);
    }
  }
  for(auto it: testVec) {
    neededCuts.loadCuts(it->info->ePos);
  }
  
  if(!isData and distats["Run"].bfind("ApplyZBoostSF") and isVSample){
    neededCuts.loadCuts(CUTS::eGen);
    neededCuts.loadCuts(CUTS::eGZ);
    neededCuts.loadCuts(CUTS::eGW);
  }
  
  neededCuts.loadCuts(_Jet->findExtraCuts());
  if(doSystematics) {
    neededCuts.loadCuts(CUTS::eGen);
  }

  for(auto it: jetCuts) {
    if(!neededCuts.isPresent(it)) continue;
    neededCuts.loadCuts(_Jet->overlapCuts(it));
  }

  if(neededCuts.isPresent(CUTS::eRWjet)) {
    neededCuts.loadCuts(_FatJet->findExtraCuts());
    neededCuts.loadCuts(_FatJet->overlapCuts(CUTS::eRWjet));
  } else {
    cout<<"WJets not needed. They will be deactivated!"<<endl;
    _FatJet->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRTau1) || neededCuts.isPresent(CUTS::eRTau2) ) {
    neededCuts.loadCuts(_Tau->findExtraCuts());
  } else {
    cout<<"Taus not needed. They will be deactivated!"<<endl;
    _Tau->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRElec1) || neededCuts.isPresent(CUTS::eRElec2) ) {
    neededCuts.loadCuts(_Electron->findExtraCuts());
  } else {
    cout<<"Electrons not needed. They will be deactivated!"<<endl;
    _Electron->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRMuon1) || neededCuts.isPresent(CUTS::eRMuon2) ) {
    neededCuts.loadCuts(_Muon->findExtraCuts());
  } else {
    cout<<"Muons not needed. They will be deactivated!"<<endl;
    _Muon->unBranch();
  }

  if( !neededCuts.isPresent(CUTS::eGen) and !isData) {
    cout<<"Gen not needed. They will be deactivated!"<<endl;
    _Gen->unBranch();
 
  }
  
  cout << "Cuts being filled: " << endl;
  for(auto cut : neededCuts.getCuts()) {
    cout << enumNames.at(static_cast<CUTS>(cut)) << "   ";
  }
  cout << endl;
}


///Smears lepton only if specified and not a data file.  Otherwise, just filles up lorentz vectors
//of the data into the vector container smearP with is in each lepton object.
void Analyzer::smearLepton(Lepton& lep, CUTS eGenPos, const PartStats& stats, const PartStats& syst_stats, int syst) {
  if( isData) {
    lep.setOrigReco();
    return;
  }

  string systname = syst_names.at(syst);
  if(!lep.needSyst(syst)) return;

  if(systname=="orig" && !stats.bfind("SmearTheParticle")){
    lep.setOrigReco();
  } else {
    systematics.loadScaleRes(stats, syst_stats, systname);
    for(size_t i = 0; i < lep.size(); i++) {
      TLorentzVector lepReco = lep.RecoP4(i);
      TLorentzVector genVec =  matchLeptonToGen(lepReco, lep.pstats["Smear"],eGenPos);
      systematics.shiftLepton(lep, lepReco, genVec, _MET->systdeltaMEx[syst], _MET->systdeltaMEy[syst], syst);
    }
  }
}

///Same as smearlepton, just jet specific
void Analyzer::smearJet(Particle& jet, const CUTS eGenPos, const PartStats& stats, int syst) {
  //at the moment
  if(isData || jet.type != PType::Jet || !stats.bfind("SmearTheJet")) {
    jet.setOrigReco();
    return;
  }
  //add energy scale uncertainty

  string systname = syst_names.at(syst);
  //cout << systname << endl;
  for(size_t i=0; i< jet.size(); i++) {
    TLorentzVector jetReco = jet.RecoP4(i);

    if(JetMatchesLepton(*_Muon, jetReco, stats.dmap.at("MuonMatchingDeltaR"), CUTS::eGMuon) ||
       JetMatchesLepton(*_Tau, jetReco, stats.dmap.at("TauMatchingDeltaR"), CUTS::eGTau) ||
       JetMatchesLepton(*_Electron, jetReco,stats.dmap.at("ElectronMatchingDeltaR"), CUTS::eGElec)){
      jet.addP4Syst(jetReco,syst);
      continue;
    }

    double sf=1.;
    //only apply corrections for jets not for FatJets

    TLorentzVector genJet=matchJetToGen(jetReco, jet.pstats["Smear"],eGenPos);
    //cout << "1: " << genJet.Pt() << endl; //01.05.18
    if(systname=="orig"){
      sf= 1.0 ; //jetScaleRes.GetRes(jetReco,genJet, rho, 0); 01.15.18
    }else if(systname=="Jet_Res_Up"){
      sf=jetScaleRes.GetRes(jetReco,genJet, rho, 1);
    }else if(systname=="Jet_Res_Down"){
      sf=jetScaleRes.GetRes(jetReco,genJet, rho, -1);
    }else if(systname=="Jet_Scale_Up"){
      sf = 1.05 ; //1.+ jetScaleRes.GetScale(jetReco, false, +1.); 01.15.18
    }else if(systname=="Jet_Scale_Down"){
      sf = 0.95 ; //1.- jetScaleRes.GetScale(jetReco, false, -1) ; 01.15.18
    }

    systematics.shiftParticle(jet, jetReco, sf, _MET->systdeltaMEx[syst], _MET->systdeltaMEy[syst], syst); 

  }
}


/////checks if jet is close to a lepton and the lepton is a gen particle, then the jet is a lepton object, so
//this jet isn't smeared
bool Analyzer::JetMatchesLepton(const Lepton& lepton, const TLorentzVector& jetV, double partDeltaR, CUTS eGenPos) {
  for(size_t j = 0; j < lepton.size(); j++) {
    if(jetV.DeltaR(lepton.RecoP4(j)) < partDeltaR && matchLeptonToGen(lepton.RecoP4(j), lepton.pstats.at("Smear"), eGenPos) != TLorentzVector(0,0,0,0)) return true;
  }
  return false;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchLeptonToGen(const TLorentzVector& lvec, const PartStats& stats, CUTS ePos) {
  if(ePos == CUTS::eGTau) {
    return matchTauToGen(lvec, stats.dmap.at("GenMatchingDeltaR"));
  }
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(_Gen->p4(it)) <= stats.dmap.at("GenMatchingDeltaR")) {
      if(stats.bfind("UseMotherID") && abs(_Gen->motherpdg_id->at(it)) != stats.dmap.at("MotherID")) continue;
      return _Gen->p4(it);
    }
  }
  return TLorentzVector(0,0,0,0);
}


///Tau specific matching fucntion.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically
TLorentzVector Analyzer::matchTauToGen(const TLorentzVector& lvec, double lDeltaR) {
  TLorentzVector genVec(0,0,0,0);
  int i = 0;
  for(vec_iter it=active_part->at(CUTS::eGTau)->begin(); it !=active_part->at(CUTS::eGTau)->end();it++, i++) {
    int nu = active_part->at(CUTS::eNuTau)->at(i);
    if(nu == -1) continue;

    genVec = _Gen->p4(*it) - _Gen->p4(nu);
    if(lvec.DeltaR(genVec) <= lDeltaR) {
      return genVec;
    }
  }
  return genVec;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchJetToGen(const TLorentzVector& lvec, const PartStats& stats, CUTS ePos) {
  //for the future store gen jets
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(_Gen->p4(it)) <= stats.dmap.at("GenMatchingDeltaR")) {
      //nothing more than b quark or gluon
      if( !(abs(_Gen->pdg_id->at(it))<5 || _Gen->pdg_id->at(it)==9 ||  _Gen->pdg_id->at(it)==21) ) continue;
      return _Gen->p4(it);
    }
  }
  return TLorentzVector(0,0,0,0);
}



////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
int Analyzer::matchToGenPdg(const TLorentzVector& lvec, double minDR) {
  double _minDR=minDR;
  int found=-1;
  for(size_t i=0; i< _Gen->size(); i++) {

    if(lvec.DeltaR(_Gen->p4(i)) <=_minDR) {
      //only hard interaction
      if( _Gen->status->at(i)<10){
        found=i;
        _minDR=lvec.DeltaR(_Gen->p4(i));
      }
    }
  }
  if (found>=0){
    return _Gen->pdg_id->at(found);
  }
  return 0;
}


////Calculates the number of gen particles.  Based on id number and status of each particle
void Analyzer::getGoodGen(const PartStats& stats) {
  if(! neededCuts.isPresent(CUTS::eGen)) return;
  for(size_t j = 0; j < _Gen->size(); j++) {
    //we are not interested in pythia info here!
    //if(_Gen->status->at(j)>10){
      //continue;
    //}
    int id = abs(_Gen->pdg_id->at(j));
    if(genMaper.find(id) != genMaper.end() && _Gen->status->at(j) == genMaper.at(id)->status) {
      if(id == 15 && (_Gen->pt(j) < stats.pmap.at("TauPtCut").first || _Gen->pt(j) > stats.pmap.at("TauPtCut").second || abs(_Gen->eta(j)) > stats.dmap.at("TauEtaCut"))) continue;
      if(id == 15 && abs(_Gen->motherpdg_id->at(j)) != 24) continue;
      if(id == 24 && (_Gen->status->at(j) != 2 && _Gen->status->at(j) != 62)) continue;
      active_part->at(genMaper.at(id)->ePos)->push_back(j);
    }
    //something special for jet
    if( (id<5 || id==9 ||  id==21) && genMaper.find(id) != genMaper.end() && _Gen->status->at(j) == genMaper.at(5)->status) {
      active_part->at(genMaper.at(5)->ePos)->push_back(j);
    }
  }
}

////Tau neutrino specific function used for calculating the number of hadronic taus
void Analyzer::getGoodTauNu() {
  for(auto it : *active_part->at(CUTS::eGTau)) {
    bool leptonDecay = false;
    int nu = -1;
    for(size_t j = 0; j < _Gen->size(); j++) {
      if(abs(_Gen->BmotherIndex->at(j)) == (it)) {
        if( (abs(_Gen->pdg_id->at(j)) == 16) && (abs(_Gen->motherpdg_id->at(j)) == 15) && (_Gen->status->at(_Gen->BmotherIndex->at(j)) == 2) ) nu = j;
        else if( (abs(_Gen->pdg_id->at(j)) == 12) || (abs(_Gen->pdg_id->at(j)) == 14) ) leptonDecay = true;
      }
    }
    nu = (leptonDecay) ? -1 : nu;
    active_part->at(CUTS::eNuTau)->push_back(nu);
  }

}

///Function used to find the number of reco leptons that pass the various cuts.
///Divided into if blocks for the different lepton requirements.
void Analyzer::getGoodRecoLeptons(const Lepton& lep, const CUTS ePos, const CUTS eGenPos, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePos)) return;

  string systname = syst_names.at(syst);
  if(!lep.needSyst(syst)) {
    active_part->at(ePos) = goodParts[ePos];
    return;
  }

  int i = 0;
  for(auto lvec: lep) {
    bool passCuts = true;
    if (fabs(lvec.Eta()) > stats.dmap.at("EtaCut")) passCuts = passCuts && false;
    else if (lvec.Pt() < stats.pmap.at("PtCut").first || lvec.Pt() > stats.pmap.at("PtCut").second) passCuts = passCuts && false;
    
    if((lep.pstats.at("Smear").bfind("MatchToGen")) && (!isData)) {   /////check
      if(matchLeptonToGen(lvec, lep.pstats.at("Smear") ,eGenPos) == TLorentzVector(0,0,0,0)) continue;
    }
    
    for( auto cut: stats.bset) {
      if(!passCuts) break;
      else if(cut == "DoDiscrByIsolation") {
        double firstIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").first : ival(ePos) - ival(CUTS::eRTau1) + 1;
        double secondIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").second : 0;
        passCuts = passCuts && lep.get_Iso(i, firstIso, secondIso);
      }
      else if(cut == "DiscrIfIsZdecay" && lep.type != PType::Tau ) passCuts = passCuts && isZdecay(lvec, lep);
      else if(cut == "DiscrByMetDphi") passCuts = passCuts && passCutRange(absnormPhi(lvec.Phi() - _MET->phi()), stats.pmap.at("MetDphiCut"));
      else if(cut == "DiscrByMetMt") passCuts = passCuts && passCutRange(calculateLeptonMetMt(lvec), stats.pmap.at("MetMtCut"));
      /////muon cuts
      else if(lep.type == PType::Muon){
        if(cut == "DoDiscrByTightID") passCuts = passCuts && _Muon->tight->at(i);
        else if(cut == "DoDiscrBySoftID") passCuts = passCuts && _Muon->soft->at(i);
      }
      ////electron cuts
      else if(lep.type == PType::Electron){
        if(cut == "DoDiscrByVetoID") passCuts = passCuts && _Electron->isPassVeto->at(i);
        else if(cut == "DoDiscrByLooseID" ) passCuts = passCuts && _Electron->isPassLoose->at(i);
        else if(cut == "DoDiscrByMediumID") passCuts = passCuts && _Electron->isPassMedium->at(i);
        else if(cut == "DoDiscrByTightID" ) passCuts = passCuts && _Electron->isPassTight->at(i);
        else if(cut == "DoDiscrByHEEPID"  ) passCuts = passCuts && _Electron->isPassHEEPId->at(i);
      }
      else if(lep.type == PType::Tau){
        if(cut == "DoDiscrByCrackCut") passCuts = passCuts && !isInTheCracks(lvec.Eta());
        /////tau cuts
        else if(cut == "DoDiscrByLeadTrack") passCuts = passCuts && (_Tau->leadChargedCandPt->at(i) >= stats.dmap.at("LeadTrackThreshold"));
             // ----Electron and Muon vetos
        else if (cut == "DoDiscrAgainstElectron") passCuts = passCuts && _Tau->pass_against_Elec(ePos, i);
        else if (cut == "SelectTausThatAreElectrons") passCuts = passCuts && !_Tau->pass_against_Elec(ePos, i);

        else if (cut == "DoDiscrAgainstMuon") passCuts = passCuts && _Tau->pass_against_Muon(ePos, i);
        else if (cut == "SelectTausThatAreMuons") passCuts = passCuts &&  !_Tau->pass_against_Muon(ePos, i);

        else if(cut == "DiscrByProngType") {
          passCuts = passCuts && (stats.smap.at("ProngType").find("hps") == string::npos || _Tau->decayModeFindingNewDMs->at(i) != 0);
          passCuts = passCuts && passProng(stats.smap.at("ProngType"), _Tau->nProngs->at(i));
        }
              // ----anti-overlap requirements
        else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      }
      else cout << "cut: " << cut << " not listed" << endl;
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }

  return;
}

////Jet specific function for finding the number of jets that pass the cuts.
//used to find the nubmer of good jet1, jet2, central jet, 1st and 2nd leading jets and bjet.
void Analyzer::getGoodRecoJets(CUTS ePos, const PartStats& stats, const int syst) {

  if(! neededCuts.isPresent(ePos)) return;

  string systname = syst_names.at(syst);
  if(!_Jet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;

  for(auto lvec: *_Jet) {
    bool passCuts = true;
    if( ePos == CUTS::eRCenJet) passCuts = passCuts && (fabs(lvec.Eta()) < 2.5);
    else  passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats.pmap.at("EtaCut"));
    passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut")) ;

    for( auto cut: stats.bset) {
      if(!passCuts) break;

    /// BJet specific
      else if(cut == "ApplyJetBTagging") passCuts = passCuts && (_Jet->bDiscriminator->at(i) > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "MatchBToGen") passCuts = passCuts && (isData ||  abs(_Jet->partonFlavour->at(i)) == 5);
      else if(cut == "ApplyLooseID") passCuts = passCuts && _Jet->passedLooseJetID(i);

    // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));

      else if(cut == "UseBtagSF") {
        double bjet_SF = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, lvec.Eta(), lvec.Pt());
        passCuts = passCuts && (isData || ((double) rand()/(RAND_MAX)) <  bjet_SF);
      }
    }
    if(_Jet->pstats["BJet"].bfind("RemoveBJetsFromJets") and ePos!=CUTS::eRBJet){
      passCuts = passCuts && find(active_part->at(CUTS::eRBJet)->begin(), active_part->at(CUTS::eRBJet)->end(), i) == active_part->at(CUTS::eRBJet)->end();
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }

  //clean up for first and second jet
  //note the leading jet has to be selected fist!
  if(ePos == CUTS::eR1stJet || ePos == CUTS::eR2ndJet) {
    vector<pair<double, int> > ptIndexVector;
    for(auto it : *active_part->at(ePos)) {
      ptIndexVector.push_back(make_pair(_Jet->pt(it),it));
    }
    sort(ptIndexVector.begin(),ptIndexVector.end());
    if(ePos == CUTS::eR1stJet && ptIndexVector.size()>0)
      active_part->at(ePos)->push_back(ptIndexVector.back().second);
    if(ePos == CUTS::eR2ndJet && ptIndexVector.size()>1)
      active_part->at(ePos)->push_back(ptIndexVector.at(ptIndexVector.size()-2).second);
  }
}


////FatJet specific function for finding the number of V-jets that pass the cuts.
void Analyzer::getGoodRecoFatJets(CUTS ePos, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePos)) return;

  string systname = syst_names.at(syst);
  if(!_FatJet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;

  for(auto lvec: *_FatJet) {
    bool passCuts = true;
    passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats.pmap.at("EtaCut"));
    passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut")) ;

    ///if else loop for central jet requirements
    for( auto cut: stats.bset) {
      if(!passCuts) break;

      else if(cut == "ApplyJetWTagging") passCuts = passCuts && (passCutRange(_FatJet->tau2->at(i)/_FatJet->tau1->at(i), stats.pmap.at("JetTau2Tau1Ratio")) &&
                                                     passCutRange(_FatJet->PrunedMass->at(i), stats.pmap.at("JetWmassCut")));
    // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));
    
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }
}

///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(overlapper.p4(it)) < MatchingDeltaR) return true;
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


///sees if the event passed one of the two cuts provided
void Analyzer::TriggerCuts(vector<int>& prevTrig, const vector<string>& trigvec, CUTS ePos) {
  if(! neededCuts.isPresent(ePos)) return;
  //cout<<" trigger "<<Trigger_decision->size()<<endl;
  if(version==1){
    for(size_t i = 0; i < trigvec.size(); i++) {
      for(size_t j =0; j<Trigger_decisionV1->size();  j++){
        //cout<<"i:  "<<prevTrig.at(i)<<" j:  "<<j<<" dec(j):  "<<Trigger_decisionV1->at(j)<<endl;
        if(prevTrig.at(i)==Trigger_decisionV1->at(j)){
          active_part->at(ePos)->push_back(0);
          return;
        }
      }
    }
  }else{
    for(int i = 0; i < (int)trigvec.size(); i++) {
      if(Trigger_decision->at(prevTrig.at(i)) == 1) {
        active_part->at(ePos)->push_back(0);
        return;
      }
    }
  }
}


////VBF specific cuts dealing with the leading jets.
void Analyzer::VBFTopologyCut(const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(CUTS::eSusyCom)) return;
  string systname = syst_names.at(syst);

  
  if(systname!="orig"){
    //only jet stuff is affected
    //save time to not rerun stuff
    if( systname.find("Jet")==string::npos){
      active_part->at(CUTS::eSusyCom)=goodParts[CUTS::eSusyCom];
      return;
    }
  }

  if(active_part->at(CUTS::eR1stJet)->size()==0 || active_part->at(CUTS::eR2ndJet)->size()==0) return;

  TLorentzVector ljet1 = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0));
  TLorentzVector ljet2 = _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));
  TLorentzVector dijet = ljet1 + ljet2;
  double dphi1 = normPhi(ljet1.Phi() - _MET->phi());
  double dphi2 = normPhi(ljet2.Phi() - _MET->phi());
  
  bool passCuts = true;
  for(auto cut: stats.bset) {
    if(!passCuts) break;
    else if(cut == "DiscrByMass") passCuts = passCuts && passCutRange(dijet.M(), stats.pmap.at("MassCut"));
    else if(cut == "DiscrByPt") passCuts = passCuts && passCutRange(dijet.Pt(), stats.pmap.at("PtCut"));
    else if(cut == "DiscrByDeltaEta") passCuts = passCuts && passCutRange(abs(ljet1.Eta() - ljet2.Eta()), stats.pmap.at("DeltaEtaCut"));
    else if(cut == "DiscrByDeltaPhi") passCuts = passCuts && passCutRange(absnormPhi(ljet1.Phi() - ljet2.Phi()), stats.pmap.at("DeltaPhiCut"));  
    else if(cut == "DiscrByOSEta") passCuts = passCuts && (ljet1.Eta() * ljet2.Eta() < 0);
    else if(cut == "DiscrByR1") passCuts = passCuts && passCutRange(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0)), stats.pmap.at("R1Cut"));
    else if(cut == "DiscrByR2") passCuts = passCuts && passCutRange(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), stats.pmap.at("R2Cut"));
    else if(cut == "DiscrByAlpha") {
      double alpha = (dijet.M() > 0) ? ljet2.Pt() / dijet.M() : -1;
      passCuts = passCuts && passCutRange(alpha, stats.pmap.at("AlphaCut"));
    }
    else if(cut == "DiscrByDphi1") passCuts = passCuts && passCutRange(abs(dphi1), stats.pmap.at("Dphi1Cut"));
    else if(cut == "DiscrByDphi2") passCuts = passCuts && passCutRange(abs(dphi2), stats.pmap.at("Dphi2Cut"));

    else cout << "cut: " << cut << " not listed" << endl;
  }

  if(passCuts)  active_part->at(CUTS::eSusyCom)->push_back(0);
  return;
}

bool Analyzer::passCutRange(double value, const pair<double, double>& cuts) {
  return (value > cuts.first && value < cuts.second);
}

//-----Calculate lepton+met transverse mass
double Analyzer::calculateLeptonMetMt(const TLorentzVector& Tobj) {
  double px = Tobj.Px() + _MET->px();
  double py = Tobj.Py() + _MET->py();
  double et = Tobj.Et() + _MET->energy();
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

  if(howCalc == "InvariantMass") {
    return (Tobj1 + Tobj2).M();
  }

  
  //////check this equation/////
  if(howCalc == "CollinearApprox") {
    double denominator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1 = (Tobj2.Py()*_MET->px() - Tobj2.Px()*_MET->py())/denominator;
    double x2 = (Tobj1.Px()*_MET->py() - Tobj1.Py()*_MET->px())/denominator;
    ratioNotInRange=!((x1 < 0.) && (x2 < 0.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px()*(1 + x1) + Tobj2.Px()*(1+x2)), (Tobj1.Py()*(1+x1) + Tobj2.Py()*(1+x2)), (Tobj1.Pz()*(1+x1) + Tobj2.Pz()*(1+x2)), (Tobj1.Energy()*(1+x1) + Tobj2.Energy()*(1+x2)) );
      return The_LorentzVect.M();
    }
  }

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    return (Tobj1 + Tobj2 + _MET->p4()).M();
  }

  return (Tobj1 + Tobj2).M();
}

////Tests if the CollinearApproximation works for finding the mass of teh particles
bool Analyzer::passDiParticleApprox(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string howCalc) {
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + _MET->px())) - (Tobj2.Px() * (Tobj1.Py() + _MET->py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + _MET->py())) - (Tobj1.Py() * (Tobj2.Px() + _MET->px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    return (x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.);
  } else {
    return true;
  }
}


/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonCombos(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePosFin)) return;
  string systname = syst_names.at(syst);

  if(!lep1.needSyst(syst) && !lep2.needSyst(syst)) {
    active_part->at(ePosFin)=goodParts[ePosFin];
    return;
  }

  bool sameParticle = (&lep1 == &lep2);
  TLorentzVector part1, part2;

  for(auto i1 : *active_part->at(ePos1)) {
    part1 = lep1.p4(i1);
    for(auto i2 : *active_part->at(ePos2)) {
      if(sameParticle && i2 <= i1) continue;
      part2 = lep2.p4(i2);
      bool passCuts = true;
      for(auto cut : stats.bset) {
        if(!passCuts) break;
        else if (cut == "DiscrByDeltaR") passCuts = passCuts && (part1.DeltaR(part2) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && passCutRange(cos(absnormPhi(part1.Phi() - part2.Phi())), stats.pmap.at("CosDphiCut"));
        else if(cut == "DiscrByDeltaPt") passCuts = passCuts && passCutRange(part1.Pt() - part2.Pt(), stats.pmap.at("DeltaPtCutValue"));
        else if(cut == "DiscrByCDFzeta2D") {
          double CDFzeta = stats.dmap.at("PZetaCutCoefficient") * getPZeta(part1, part2).first
            + stats.dmap.at("PZetaVisCutCoefficient") * getPZeta(part1, part2).second;
          passCuts = passCuts && passCutRange(CDFzeta, stats.pmap.at("CDFzeta2DCutValue"));
        }
        else if(cut == "DiscrByDeltaPtDivSumPt") {
          double ptDiv = (part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt());
          passCuts = passCuts && passCutRange(ptDiv, stats.pmap.at("DeltaPtDivSumPtCutValue"));
        }
        else if (cut == "DiscrByMassReco") {
          double diMass = diParticleMass(part1,part2, stats.smap.at("HowCalculateMassReco"));
          passCuts = passCuts && passCutRange(diMass, stats.pmap.at("MassCut"));
        }
        else if(cut == "DiscrByCosDphiPtAndMet"){       
          double CosDPhi1 = cos(absnormPhi(part1.Phi() - _MET->phi()));       
          passCuts = passCuts && passCutRange(CosDPhi1, stats.pmap.at("CosDphiPtAndMetCut"));
        }


        else cout << "cut: " << cut << " not listed" << endl;
      }

       if (stats.bfind("DiscrByOSLSType")){
          //   if it is 1 or 0 it will end up in the bool map!!
         if(stats.bfind("DiscrByOSLSType") && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
       }else if (stats.dmap.find("DiscrByOSLSType") != stats.dmap.end() ){
         if(lep1.charge(i1) * lep2.charge(i2) > 0) continue;
       }else if (stats.smap.find("DiscrByOSLSType") != stats.smap.end() ){
         if(stats.smap.at("DiscrByOSLSType") == "LS" && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
         else if(stats.smap.at("DiscrByOSLSType") == "OS" && (lep1.charge(i1) * lep2.charge(i2) >= 0)) continue;
       }

      ///Particlesp that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2
      if(passCuts) active_part->at(ePosFin)->push_back(i1*BIG_NUM + i2);
    }
  }
}


//////////////LOOK INTO DIJET PICKING
///////HOW TO GET RID OF REDUNCENCIES??

/////Same as gooddilepton, just jet specific
void Analyzer::getGoodDiJets(const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(CUTS::eDiJet)) return;
  string systname = syst_names.at(syst);
  if(systname!="orig"){
    //save time to not rerun stuff
    if( systname.find("Jet")==string::npos){
      active_part->at(CUTS::eDiJet)=goodParts[CUTS::eDiJet];
      return;
    }
  }
  TLorentzVector jet1, jet2;
  // ----Separation cut between jets (remove overlaps)
  for(auto ij2 : *active_part->at(CUTS::eRJet2)) {
    jet2 = _Jet->p4(ij2);
    for(auto ij1 : *active_part->at(CUTS::eRJet1)) {
      if(ij1 == ij2) continue;
      jet1 = _Jet->p4(ij1);
      
      bool passCuts = true;
      //cout<<"---------------------"<<endl;
      for(auto cut : stats.bset) {
        //cout<<cut<<"    "<<passCuts<<endl;;
        if(!passCuts) break;
        else if(cut == "DiscrByDeltaR") passCuts = passCuts && (jet1.DeltaR(jet2) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByDeltaEta") passCuts = passCuts && passCutRange(abs(jet1.Eta() - jet2.Eta()), stats.pmap.at("DeltaEtaCut"));
        else if(cut == "DiscrByDeltaPhi") passCuts = passCuts && passCutRange(absnormPhi(jet1.Phi() - jet2.Phi()), stats.pmap.at("DeltaPhiCut"));  
        else if(cut == "DiscrByOSEta") passCuts = passCuts && (jet1.Eta() * jet2.Eta() < 0);
        else if(cut == "DiscrByMassReco") passCuts = passCuts && passCutRange((jet1+jet2).M(), stats.pmap.at("MassCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && passCutRange(cos(absnormPhi(jet1.Phi() - jet2.Phi())), stats.pmap.at("CosDphiCut"));
        else cout << "cut: " << cut << " not listed" << endl;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2
      if(passCuts) active_part->at(CUTS::eDiJet)->push_back(ij1*_Jet->size() + ij2);
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
  for(vector<TLorentzVector>::const_iterator lepit= lep.begin(); lepit != lep.end(); lepit++) {
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
  double px = visPx + _MET->px();
  double py = visPy + _MET->py();
  return make_pair(px*zetaX + py*zetaY, visPx*zetaX + visPy*zetaY);
}

//---MY WEIGHTING and DECAY LIST----------------------------------------------------------------------------------------------------------------------------------------
// double Analyzer::getZBoostWeight(){  //new7.28.17
//   double z2_p = 0.; //new9.14.17
//   for (unsigned p=0; p < _Gen->pt->size(); p++){ //new9.14.17 Loop over the full size of _Gen
//     if (_Gen->pdg_id->at(p) == 23 && (_Gen->status->at(p) == 2 || _Gen->status->at(p) == 62)){ //new9.14.17 If we get a Z with the right status...
//       z2_p = _Gen->pt->at(p); //new9.14.17 ...we grab its p_T...
//       break; //new9.14.17 ..and then get out of the loop once we have it.
//     }
//   }
//   return z2_p;
//   //return (Tobj1 + Tobj2).Pt();  //new7.28.17
// }//new7.28.17

// double Analyzer::getWBoostWeight(){  //new7.28.17
//   double w2_p = 0.;  //new7.28.17
//   for (unsigned p=0; p < _Gen->pt->size(); p++){  //new7.28.17 Loop over the full size of _Gen.
//     if (_Gen->pdg_id->at(p) == 24 && (_Gen->status->at(p) == 2 || _Gen->status->at(p) == 62)){  //new7.28.17 If we get a W with the right status...
//     w2_p = _Gen->pt->at(p);  //new7.28.17 ...we grab its p_T...
//     break;  //new7.28.17 ...and then get out of the loop once we have it.
//     }  //new7.28.17
//   }  //new7.28.17
//     return w2_p;  //new7.28.17 We want to return the p_T of the W.
// }  //new7.28.17
  
// double Analyzer::getParallelRecoilWeight(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2){
//   //-----Get the recoil Lorentz vector.-----//
//   double boostzpx = Tobj1.Px() + Tobj2.Px();     //x-component of zboost
//   double boostzpy = Tobj1.Py() + Tobj2.Py();     //y-component of zboost
//   TLorentzVector boostzplorentz;     //Make Lorentz vector for zboost.           
//   boostzplorentz.SetPxPyPzE(boostzpx, boostzpy, 0, 0);     //Fill Lorentz vec. 
//   TLorentzVector u_t;     //Make Lorentz vector for recoil.
//   u_t = -1 * (boostzplorentz + theMETVector);     //Fill Lorentz vec (will use Px, Py).
    
//   //-----2-vector for zboost now.-----//
//   TVector2 boostzp;     //two-component boost (Pt) vector for easy manipulations
//   boostzp.Set(boostzpx, boostzpy);     //Fill TVector 2.

//   //-----Get unit vectors parallel and perpendicular to boost.-----//
//   TVector2 boostzp_unitpar;     //two-component unit vector parallel to the Z Pt
//   boostzp_unitpar = boostzp.Unit();
//   TVector2 boostzp_unitper;     //two-component unit vector perp. to the Z Pt
//   boostzp_unitper.Set((-1 / boostzp.Mod()) * boostzpy, (1 / boostzp.Mod()) * boostzpx);

//   //-----Project the recoil onto the parallel unit vector to zboost.-----//
//   TVector2 ut_pt2;     //two-component vector for the recoil Pt in each direction.
//   ut_pt2.Set(u_t.Px(), u_t.Py());
//   double ut_par = boostzp_unitpar * ut_pt2;     //Project recoil onto unit vector parallel to the Z Pt.
  
//   return ut_par;
// }

// void Analyzer::checkParticleDecayList(int event){ //new7.25.17
//   BOOM->GetEntry(event); //new7.25.17
//   fstream file;  //new7.25.17
//   file.open("particle_decay_list.txt", ios_base::out | ios_base::in);  //will not create file //new7.25.17
//   int s;  //new7.25.17
//   if (file.is_open()){ //new7.25.17
//     cout << "Warning, file already exists.";  //new7.25.17
//     cout << "Do you wish to clear the file? 1 for yes; 0 for no.";  //new7.25.17
//     cin >> s; //new7.25.17
//     if (s == 1)  //new7.25.17
//       {file.open("particle_decay_list.txt", ios::out | ios::trunc); //new7.25.17
// 	file.close(); //new7.25.17
// 	cout << "You have cleared the file.";} //new7.25.17
//   }
// }


// void Analyzer::writeParticleDecayList(int event, fstream& file){  //new7.25.17
//   BOOM->GetEntry(event); //new7.25.17
//   bool isThereATau = false;  //new7.25.17
//   for (unsigned s=0; s < _Gen->pt->size(); s++){  //new7.25.17
//     if (_Gen->pdg_id->at(s) == 15){  //new7.25.17
//       isThereATau = true;  //new7.25.17
//       break;  //new7.25.17
//     }  //new7.25.17
//   }  //new7.25.17
//   if(!isThereATau)  return;  //new7.25.17
//   for (unsigned p=0; p < _Gen->pt->size(); p++){  //new7.25.17
//     file << setw(2) << p << setw(2) << " " << setw(8) << "pdg_id: " << setw(4) << abs(_Gen->pdg_id->at(p)) << setw(2) << "  " << setw(5) << "p_T: " << setw(10) << _Gen->pt->at(p) << setw(2) << "  " << setw(5) << "phi: " << setw(10) << _Gen->phi->at(p) << setw(10) << "status: " << setw(3) << _Gen->status->at(p) << setw(2) << "  " << setw(7) << "mind:" << setw(1) << " " << setw(2) << _Gen->BmotherIndex->at(p) << "\n";  //new7.25.17
//   } //new7.25.17
//   return;  //new7.25.17
// } //new7.25.17

// bool Analyzer::checkforW(){  //new7.28.17 Function checks to see if we have the W we want in _Gen.
//   bool isThereAW = false;  //new7.28.17 Initiate the boolean "isThereAW" to false.  

//   for (unsigned s=0.; s < _Gen->pt->size(); s++){  //new7.28.17 Over the full size of _Gen...
//     if (abs(_Gen->pdg_id->at(s)) == 24 && (_Gen->status->at(s) == 62 || _Gen->status->at(s) == 2) ){  //new7.28.17 If the particle is a W with the appropriate status...
//       isThereAW = true;  //new7.28.17 ...toggle the boolean to true.
//       break;  //new7.28.17 Get out of the loop over _Gen once you get a W.
//     }
//   }

//   return isThereAW;  //new7.28.17 Get out whether we have the W with this bool.
// }

// bool Analyzer::checkforZ(){  //new9.14.17 Function checks to see if we have the Z we want in _Gen.
//   bool isThereAZ = false;  //new9.14.17 Initiate the boolean "isThereAZ" to false.                                                                                    
//   for (unsigned s=0.; s < _Gen->pt->size(); s++){  //new9.14.17 Over the full size of _Gen...                                                                                
//     if (abs(_Gen->pdg_id->at(s)) == 23 && (_Gen->status->at(s) == 62 || _Gen->status->at(s) == 2) ){  //new9.14.17 If the particle is a Z with the appropriate status...
//       isThereAZ = true;  //new9.14.17 ...toggle the boolean to true.
//       break;  //new9.14.17 Get out of the loop over _Gen once you get a Z.                                                                                              
//     }
//   }

//   return isThereAZ;  //new9.14.17 Get out whether we have the W with this bool.                                                                                          
// }
//---END MY WEIGHTING and DECAY LIST------------------------------------------------------------------------------------------------------------------------------------

double Analyzer::getZBoostWeight(){
  double boostweigth=1.;
  if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGZ)->size() ==1 || active_part->at(CUTS::eGW)->size() ==1)){
    //cout<<" Z or W " <<endl;
    double boostz = 0;
    if(active_part->at(CUTS::eGZ)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGZ)->at(0));
    }
    if(active_part->at(CUTS::eGW)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGW)->at(0));
    }
    if(boostz > 0 && boostz <= 50) {boostweigth = 1.1192;}
    else if (boostz > 50 && boostz <= 100) {boostweigth = 1.1034;}
    else if (boostz > 100 && boostz <= 150) {boostweigth = 1.0675;}
    else if (boostz > 150 && boostz <= 200) {boostweigth = 1.0637;}
    else if (boostz > 200 && boostz <= 300) {boostweigth = 1.0242;}
    else if (boostz > 300 && boostz <= 400) {boostweigth = 0.9453;}
    else if (boostz > 400 && boostz <= 600) {boostweigth = 0.8579;}
    else if (boostz >= 600) {boostweigth = 0.7822;}
    else {boostweigth = 1;}
  }
  return boostweigth;
}

////Grabs a list of the groups of histograms to be filled and asked Fill_folder to fill up the histograms
void Analyzer::fill_histogram() {
  if(distats["Run"].bfind("ApplyGenWeight") && gen_weight == 0.0) return;

  if(isData && blinded && maxCut == SignalRegion) return;

  const vector<string>* groups = histo.get_groups();
  if(!isData){
    wgt = 1.;
    if(distats["Run"].bfind("UsePileUpWeight")) wgt*= pu_weight;
    if(distats["Run"].bfind("ApplyGenWeight")) wgt *= (gen_weight > 0) ? 1.0 : -1.0;
    //add weight here
    if(distats["Run"].bfind("ApplyTauIDSF")) wgt *= getTauDataMCScaleFactor(0);

    if(distats["Run"].bfind("ApplyZBoostSF") && isVSample){
      wgt *= getZBoostWeight();
    }
  }else  wgt=1.;
  //backup current weight
  backup_wgt=wgt;
  cout << "syst_names size: " << syst_names.size() << endl;
  for(size_t i = 0; i < syst_names.size(); i++) { 
    for(Particle* ipart: allParticles) ipart->setCurrentP(i);
    _MET->setCurrentP(i);
    active_part =&syst_parts.at(i);
    //////i == 0 is orig or no syst case
    if(i == 0) {
      active_part = &goodParts;
      fillCuts(true);
      for(auto it: *groups) {
	cout << "it1: " << it << endl;
        fill_Folder(it, maxCut, histo, false);
      }
    }else{
      //cout << "isyst:  " << i << endl;
      //get the non particle conditions:
      for(auto itCut : nonParticleCuts){
	//cout << "itcut: " << itCut << endl;
        active_part->at(itCut)=goodParts.at(itCut);
      }
      //cout<<"________________"<<i<<endl;
      if(!fillCuts(false)) continue;
      cout << "i: " << i << endl;
      for(auto it: *syst_histo.get_groups()) { //categories in the hist info. files
	cout << "it2: " << it << endl;//cout << "isyst: " << i << endl;
        fill_Folder(it, i, syst_histo, true);
	//cout << "i2: " << i << endl;
      }
    }
  }
}


///Function that fills up the histograms
void Analyzer::fill_Folder(string group, const int max, Histogramer &ihisto, bool issyst) {
  /*be aware in this function
   * the following definition is used:
   * histAddVal(val, name) histo.addVal(val, group, max, name, wgt)
   * so each histogram knows the group, max and weight!
   */
  if(group == "FillRun" && (&ihisto==&histo)) {
    if(crbins != 1) {
      for(int i = 0; i < crbins; i++) {
        ihisto.addVal(false, group, i, "Events", 1);
        if(distats["Run"].bfind("ApplyGenWeight")) {
          //put the weighted events in bin 3
          ihisto.addVal(2, group,i, "Events", (gen_weight > 0) ? 1.0 : -1.0);
        }
        ihisto.addVal(wgt, group, i, "Weight", 1);
      }
    }
    else{
      ihisto.addVal(false, group,ihisto.get_maxfolder(), "Events", 1);
      if(distats["Run"].bfind("ApplyGenWeight")) {
        //put the weighted events in bin 3
        ihisto.addVal(2, group,ihisto.get_maxfolder(), "Events", (gen_weight > 0) ? 1.0 : -1.0);
      }
      ihisto.addVal(wgt, group, ihisto.get_maxfolder(), "Weight", 1);
    }
    histAddVal(true, "Events");
    histAddVal(bestVertices, "NVertices");
  } else if(!isData && group == "FillGen") {

    int nhadtau = 0;
    TLorentzVector genVec;
    int i = 0;
    for(vec_iter it=active_part->at(CUTS::eGTau)->begin(); it!=active_part->at(CUTS::eGTau)->end(); it++, i++) {
      int nu = active_part->at(CUTS::eNuTau)->at(i);
      if(nu != -1) {
        genVec = _Gen->p4(*it) - _Gen->p4(nu);
        histAddVal(genVec.Pt(), "HadTauPt");
        histAddVal(genVec.Eta(), "HadTauEta");
        nhadtau++;
	//--MY CODE FROM VALIDATING ALEJANDRO'S W-STUDIES---------------------------------------------------------------------------------------------------------------
	//double ghadtau_frac_scalar = Hadtau.Pt() / (Hadtau.Pt()+nu.Pt());
	//double ghadtau_frac_vector = Hadtau.Pt() / (Gentau.Pt());
	//histAddVal(ghadtau_frac_scalar, "HadTauFrac_scalar");
	//histAddVal(ghadtau_frac_vector, "HadTauFrac_vector");
	//histAddVal(Hadtau.Pt(), "HadTauPt");
	//histAddVal(nu.Pt(), "NuPt");
	//histAddVal(nu.P(), "NuP");
	//histAddVal(Hadtau.Eta(), "HadTauEta");
	//histAddVal(Hadtau.DeltaR(nu), "DeltaR");
	//---END MY CODE FROM VALIDATING ALEJANDRO'S W-STUDIES----------------------------------------------------------------------------------------------------------
      }
      histAddVal(_Gen->energy(*it), "TauEnergy");
      histAddVal(_Gen->pt(*it), "TauPt");
      histAddVal(_Gen->eta(*it), "TauEta");
      histAddVal(_Gen->phi(*it), "TauPhi");
      for(vec_iter it2=it+1; it2!=active_part->at(CUTS::eGTau)->end(); it2++) {
        histAddVal(diParticleMass(_Gen->p4(*it),_Gen->p4(*it2), "none"), "DiTauMass");
      }
    }
    histAddVal(active_part->at(CUTS::eGTau)->size(), "NTau");
    histAddVal(nhadtau, "NHadTau");

    for(auto it : *active_part->at(CUTS::eGZ)) {
      histAddVal(_Gen->pt(it), "ZPt");
      histAddVal(_Gen->eta(it), "ZEta");
      histAddVal(_Gen->p4(it).M(), "ZMass");
    }
    histAddVal(active_part->at(CUTS::eGZ)->size(), "NZ");

    for(auto it : *active_part->at(CUTS::eGW)) {
      histAddVal(_Gen->pt(it), "WPt");
      histAddVal(_Gen->eta(it), "WEta");
      histAddVal(_Gen->p4(it).M(), "WMass");
    }
    histAddVal(active_part->at(CUTS::eGW)->size(), "NW");

    for(auto it : *active_part->at(CUTS::eGMuon)) {
      histAddVal(_Gen->energy(it), "MuonEnergy");
      histAddVal(_Gen->pt(it), "MuonPt");
      histAddVal(_Gen->eta(it), "MuonEta");
      histAddVal(_Gen->phi(it), "MuonPhi");
    }
    histAddVal(active_part->at(CUTS::eGMuon)->size(), "NMuon");

    double mass=0;
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      //if a Z boson is explicitly there
      if(abs(_Gen->pdg_id->at(i))==11 or abs(_Gen->pdg_id->at(i))==13 or abs(_Gen->pdg_id->at(i))==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          mass=(lep1+lep2).M();
          //cout<<"mass  leptons "<<mass<<endl;
          break;
        }else{
          //cout<<_Gen->size()<<"   "<<igen<<endl;
          //if(_Gen->size()>_Gen->cur_P.size()){
           //_Gen->init();
          //}
          lep1= _Gen->RecoP4(i);
        }
      }
    }
    histAddVal(mass, "LeptonMass");
  } else if(fillInfo[group]->type == FILLER::Single) {
    Particle* part = fillInfo[group]->part;
    CUTS ePos = fillInfo[group]->ePos;

    for(auto it : *active_part->at(ePos)) {
      histAddVal(part->p4(it).Energy(), "Energy");
      histAddVal(part->p4(it).Pt(), "Pt");
      histAddVal(part->p4(it).Eta(), "Eta");
      histAddVal(part->p4(it).Phi(), "Phi");
      histAddVal(part->p4(it).DeltaPhi(_MET->p4()), "MetDphi");

      if(part->type == PType::Tau) {
        if(_Tau->nProngs->at(it) == 1){
          histAddVal(part->pt(it), "Pt_1prong");
        }else if(_Tau->nProngs->at(it) == 3){
          histAddVal(part->pt(it), "Pt_3prong");
        }
        histAddVal(_Tau->nProngs->at(it), "NumSignalTracks");
        histAddVal(_Tau->charge(it), "Charge");
        histAddVal(_Tau->leadChargedCandPt->at(it), "SeedTrackPt");
	//histAddVal(calculateLeptonMetMt(part->p4(it)), "MetMt");
      }
      if(part->type != PType::Jet) {
        histAddVal(calculateLeptonMetMt(part->p4(it)), "MetMt");
      }
      if(part->type == PType::FatJet ) {
        histAddVal(_FatJet->PrunedMass->at(it), "PrunedMass");
        histAddVal(_FatJet->SoftDropMass->at(it), "SoftDropMass");
        histAddVal(_FatJet->tau1->at(it), "tau1");
        histAddVal(_FatJet->tau2->at(it), "tau2");
        histAddVal(_FatJet->tau2->at(it)/_FatJet->tau1->at(it), "tau2Overtau1");
      }
    }

    if((part->type != PType::Jet ) && active_part->at(ePos)->size() > 0) {
      vector<pair<double, int> > ptIndexVector;
      for(auto it : *active_part->at(ePos)) {
        ptIndexVector.push_back(make_pair(part->pt(it),it));
      }
      sort(ptIndexVector.begin(),ptIndexVector.end());
      if(ptIndexVector.size()>0){
        histAddVal(part->pt(ptIndexVector.back().second), "FirstLeadingPt");
        histAddVal(part->eta(ptIndexVector.back().second), "FirstLeadingEta");
      }
      if(ptIndexVector.size()>1){
        histAddVal(part->pt(ptIndexVector.at(ptIndexVector.size()-2).second), "SecondLeadingPt");
        histAddVal(part->eta(ptIndexVector.at(ptIndexVector.size()-2).second), "SecondLeadingEta");
      }
    }

    histAddVal(active_part->at(ePos)->size(), "N");


  } else if(group == "FillMetCuts") {
    histAddVal(_MET->MHT(), "MHT");
    histAddVal(_MET->HT(), "HT");
    histAddVal(_MET->HT() + _MET->MHT(), "Meff");
    histAddVal(_MET->pt(), "Met");
    histAddVal(_MET->phi(), "MetPhi");

  } else if(group == "FillLeadingJet" && active_part->at(CUTS::eSusyCom)->size() == 0) {
    
    if(active_part->at(CUTS::eR1stJet)->size()>0) {
      histAddVal(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)).Pt(), "FirstPt");
      histAddVal(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)).Eta(), "FirstEta");
    }
    if(active_part->at(CUTS::eR2ndJet)->size()>0) {
      histAddVal(_Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0)).Pt(), "SecondPt");
      histAddVal(_Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0)).Eta(), "SecondEta");
    }


  } else if(group == "FillLeadingJet" && active_part->at(CUTS::eSusyCom)->size() != 0) {

    TLorentzVector first = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0));
    TLorentzVector second = _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

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
    double dphi1 = normPhi(first.Phi() - _MET->phi());
    double dphi2 = normPhi(second.Phi() - _MET->phi());
    double alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histAddVal(dphiDijets, "LeadSublDijetDphi");
    histAddVal2(_MET->pt(),dphiDijets, "MetVsDiJetDeltaPhiLeadSubl");
    histAddVal2(fabs(first.Eta()-second.Eta()), dphiDijets, "DeltaEtaVsDeltaPhiLeadSubl");

    histAddVal(absnormPhi(_MET->phi() - LeadDiJet.Phi()), "MetDeltaPhi");


    histAddVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), "R1");
    histAddVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), "R2");
    histAddVal(normPhi(first.Phi() - _MET->MHTphi()), "Dphi1MHT");
    histAddVal(normPhi(second.Phi() - _MET->MHTphi()), "Dphi2MHT");
    histAddVal(dphi1, "Dphi1");
    histAddVal(dphi2, "Dphi2");
    histAddVal2(dphi1,dphi2, "Dphi1VsDphi2");
    histAddVal(alpha, "Alpha");


  } else if(group == "FillDiJet") {
    double leaddijetmass = 0;
    double leaddijetpt = 0;
    double leaddijetdeltaR = 0;
    double leaddijetdeltaEta = 0;
    double etaproduct = 0;
    for(auto it : *active_part->at(CUTS::eDiJet)) {
      int p1 = (it) / _Jet->size();
      int p2 = (it) % _Jet->size();
      TLorentzVector jet1 = _Jet->p4(p1);
      TLorentzVector jet2 = _Jet->p4(p2);
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
    
    for(auto index : *(active_part->at(CUTS::eRTau1)) ) {
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetmass, "mTvsLeadingMass");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetdeltaEta, "mTvsLeadingDeltaEta");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetdeltaR, "mTvsLeadingDeltaR");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetpt, "mTvsLeadingPt");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetmass, "MetDphiVSLeadingMass");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetdeltaEta, "MetDphiVSLeadingDeltaEta");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetdeltaR, "MetDphiVSLeadingDeltaR");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetpt, "MetDphiVSLeadingPt");
    }



    ////diparticle stuff
  } else if(fillInfo[group]->type == FILLER::Dipart) {
    Lepton* lep1 = static_cast<Lepton*>(fillInfo[group]->part);
    Lepton* lep2 = static_cast<Lepton*>(fillInfo[group]->part2);
    CUTS ePos = fillInfo[group]->ePos;
    string digroup = group;
    digroup.erase(0,4);

    TLorentzVector part1;
    TLorentzVector part2;


    for(auto it : *active_part->at(ePos)) {

      int p1= (it) / BIG_NUM;
      int p2= (it) % BIG_NUM;

      part1 = lep1->p4(p1);
      part2 = lep2->p4(p2);
      
      double add_px = (part1.Pt() * cos(part1.Phi())) + (part2.Pt() * cos(part2.Phi()));
      double add_py = (part1.Pt() * sin(part1.Phi())) + (part2.Pt() * sin(part2.Phi()));
//    TLorentzVector add_p = part1 + part2;
      //double mag_addp2 = pow(add_px, 2) + pow(add_py, 2);
      TLorentzVector boost_pt;     //initialized vector for x and y-components of the muon momenta
      boost_pt.SetPxPyPzE(add_px, add_py, 0, 0);
      //double booster = boost_pt.Pt();                                                                                                   
      TLorentzVector u_t;     //initialized vector for recoil                                                                                                
      TLorentzVector met_pt;
      met_pt.SetPxPyPzE(_MET->px(), _MET->py(), 0, 0);
      u_t = -1 * (boost_pt + met_pt);     //Recoil is the negative sum of the Z Pt and MET.
      TVector2 boost_pt2;     //two-component boost (Pt) vector for easy manipulations                                                                                 
      boost_pt2.Set(add_px, add_py);     //put the muons' x and y-momenta in it
      TVector2 boost_unitpar;     //two-component vector for the unit vector parallel to the Z Pt
      boost_unitpar = boost_pt2.Unit();
      TVector2 boost_unitper;     //two-component vector for the unit vector perpendicular to the Z Pt                                                                  
      boost_unitper.Set((-1 / boost_pt2.Mod()) * add_py, (1 / boost_pt2.Mod()) * add_px);
      TVector2 ut_pt2;     //two-component vector for the recoil Pt                      
      ut_pt2.Set(u_t.Px(), u_t.Py());
      double ut_par = boost_unitpar * ut_pt2;     //Project the recoil onto the unit vector parallel to the Z Pt                                                       
      histAddVal(ut_par, "Utpar");     //Fill the histogram for the parallel component of the recoil.
      double ut_per = boost_unitper * ut_pt2;     //Project the recoil onto the unit vector perpendicular to the Z Pt                                   
      histAddVal(ut_per, "Utper");     //Fill the histogram for the perpendicular component of the recoil.                                                              

      //do our dirty tree stuff here: 
      //zBoostTree["treeweight"] = wgt;
      //zBoostTree["mu1pt"] = part1.Pt(); 
      //zBoostTree["mu2pt"] = part2.Pt();  
      //zBoostTree["addpxZ"] = boost_pt.Px();
      //zBoostTree["addpyZ"] = boost_pt.Py();
      //zBoostTree["utpx"] = addpxZMETx;                                                                                                                           
      //zBoostTree["utpy"] = addpyZMETy;     
      //zBoostTree["utmagn"] = utmag;                                                                                                                                   
      //zBoostTree["utpar"] = ut_par;                                                                                                                                   
      //zBoostTree["utper"] = ut_per;                                                                                                                                   
      //zBoostTree["boostedz"] = TMath::Sqrt(pow(lep1->smearP.at(p1).Px() + lep2->smearP.at(p2).Px(), 2) + pow(lep1->smearP.at(p1).Py() + lep2->smearP.at(p2).Py(), 2));
      //add this and remove the other settings                                                                                                                          
      //if(goodParts[CUTS::eRJet1]->size()>0){                                                                                                                          
      //zBoostTree["FLJetPt"] =_Jet->smearP.at(goodParts[CUTS::eRJet1]->at(0)).Pt();                                                                                    
      //}else{zBoostTree["FLJetPt"] =0;                                                                                                                                 
      //}                                                                                                                                                               
      //if(goodParts[CUTS::eDiMuon]->size()>0 && it==goodParts[ePos]->begin()){                                                                                         
      //histo.fillTree("zboost");}         

      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR");
      if(group.find("Di") != string::npos) {
        histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");  
        histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
        histAddVal(boost_pt2.Mod(), "Boost");
      } else {
        histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");  
        histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
        histAddVal(boost_pt2.Mod(), "Boost");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");

      histAddVal(cos(absnormPhi(part1.Phi() - _MET->phi())), "Part1CosDphiPtandMet");
      histAddVal(cos(absnormPhi(part2.Phi() - _MET->phi())), "Part2CosDphiPtandMet");


      histAddVal(absnormPhi(part1.Phi() - _MET->phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - _MET->phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - _MET->phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - _MET->phi())), "CosDphi_DeltaPtAndMet");


      double diMass = diParticleMass(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"));
      if(passDiParticleApprox(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"))) {
        histAddVal(diMass, "ReconstructableMass");
      } else {
        histAddVal(diMass, "NotReconstructableMass");
      }

      double InvMass = diParticleMass(part1,part2, "InvariantMass");
      histAddVal(InvMass, "InvariantMass");

      double ptSum = part1.Pt() + part2.Pt();
      histAddVal(ptSum, "SumOfPt");

      double PZeta = getPZeta(part1,part2).first;
      double PZetaVis = getPZeta(part1,part2).second;
      histAddVal(calculateLeptonMetMt(part1), "Part1MetMt");
      histAddVal(calculateLeptonMetMt(part2), "Part2MetMt");
      histAddVal(lep2->charge(p2) * lep1->charge(p1), "OSLS");
      histAddVal(PZeta, "PZeta");
      histAddVal(PZetaVis, "PZetaVis");
      histAddVal2(PZetaVis,PZeta, "Zeta2D");
      histAddVal((distats.at(digroup).dmap.at("PZetaCutCoefficient") * PZeta) + (distats.at(digroup).dmap.at("PZetaVisCutCoefficient") * PZetaVis), "Zeta1D");

      if ((active_part->at(CUTS::eR1stJet)->size()>0 && active_part->at(CUTS::eR1stJet)->at(0) != -1) && (active_part->at(CUTS::eR2ndJet)->size()>0 && active_part->at(CUTS::eR2ndJet)->at(0) != -1)) {
        TLorentzVector TheLeadDiJetVect = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)) + _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

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

void Analyzer::initializePileupInfo(string MCHisto, string DataHisto, string DataHistoName, string MCHistoName) {
  TFile *file1 = new TFile((PUSPACE+MCHisto).c_str());
  TH1D* histmc = (TH1D*)file1->FindObjectAny(MCHistoName.c_str());
  if(!histmc) throw std::runtime_error("failed to extract histogram");

  TFile* file2 = new TFile((PUSPACE+DataHisto).c_str());
  TH1D* histdata = (TH1D*)file2->FindObjectAny(DataHistoName.c_str());
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

///Normalizes phi to be between -PI and PI
double normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}


///Takes the absolute value of of normPhi (made because constant use)
double absnormPhi(double phi) {
  return abs(normPhi(phi));
}

