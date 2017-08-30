#ifndef Analyzer_h
#define Analyzer_h

struct CRTester;

// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>

#include <TDirectory.h>
#include <TEnv.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1.h>

#include "Particle.h"
#include "MET.h"
#include "Histo.h"

/////fix
#include "./btagging/BTagCalibrationStandalone.h"

#include "Cut_enum.h"
#include "FillInfo.h"
#include "CRTest.h"
#include "Systematics.h"
#include "JetScaleResolution.h"
#include "DepGraph.h"


double normPhi(double phi);
double absnormPhi(double phi);

//#define const
using namespace std;

static const int nTrigReq = 2;

class Analyzer {
  friend class CRTester;
public:
  Analyzer(vector<string>, string, bool setCR = false, string configFolder="PartDet");
  ~Analyzer();
  void clear_values();
  void preprocess(int);
  bool fillCuts(bool);
  void printCuts();
  void writeout();
  int nentries;
  void fill_histogram();
  void setControlRegions() { histo.setControlRegions();}

  vector<int>* getList(CUTS ePos) {return goodParts[ePos];}
  double getMet() {return _MET->pt();}
  double getHT() {return _MET->HT();}
  double getMHT() {return _MET->MHT();}
  double getMass(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string partName) {
    return diParticleMass(Tobj1, Tobj2, distats[partName].smap.at("HowCalculateMassReco"));
  }
  double getZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string partName) {
    return distats[partName].dmap.at("PZetaCutCoefficient") * getPZeta(Tobj1, Tobj2).first;
  }


private:
  void CRfillCuts();
  ///// Functions /////
  //void fill_Folder(string, const int, string syst="");
  void fill_Folder(string, const int, Histogramer& ihisto, bool issyst);

  void getInputs();
  void setupJob(string);
  void initializePileupInfo(string, string, string, string);
  void read_info(string);
  void setupGeneral();
  void initializeTrigger();
  void setCutNeeds();

  void smearLepton(Lepton&, CUTS, const PartStats&, const PartStats&, int syst=0);
  void smearJet(Particle&, CUTS, const PartStats&, int syst=0);

  bool JetMatchesLepton(const Lepton&, const TLorentzVector&, double, CUTS);
  TLorentzVector matchLeptonToGen(const TLorentzVector&, const PartStats&, CUTS);
  TLorentzVector matchTauToGen(const TLorentzVector&, double);
  TLorentzVector matchJetToGen(const TLorentzVector&, const PartStats&, CUTS);

  int matchToGenPdg(const TLorentzVector& lvec, double minDR);


  void getGoodParticles(int);
  void getGoodTauNu();
  void getGoodGen(const PartStats&);
  void getGoodRecoLeptons(const Lepton&, const CUTS, const CUTS, const PartStats&, const int);
  void getGoodRecoJets(CUTS, const PartStats&, const int);
  void getGoodRecoFatJets(CUTS, const PartStats&, const int);

  void getGoodLeptonCombos(Lepton&, Lepton&, CUTS,CUTS,CUTS, const PartStats&, const int);
  void getGoodDiJets(const PartStats&, const int);

  void VBFTopologyCut(const PartStats&, const int);
  void TriggerCuts(vector<int>&, const vector<string>&, CUTS);


  double calculateLeptonMetMt(const TLorentzVector&);
  double diParticleMass(const TLorentzVector&, const TLorentzVector&, string);
  bool passDiParticleApprox(const TLorentzVector&, const TLorentzVector&, string);
  bool isZdecay(const TLorentzVector&, const Lepton&);

  bool isOverlaping(const TLorentzVector&, Lepton&, CUTS, double);
  bool passProng(string, int);
  bool isInTheCracks(float);
  bool passedLooseJetID(int);
  bool select_mc_background();
  double getTauDataMCScaleFactor(int updown);


  pair<double, double> getPZeta(const TLorentzVector&, const TLorentzVector&);
  void create_fillInfo();

  inline bool passCutRange(string, double, const PartStats&);
  bool passCutRange(double, const pair<double, double>&);
  bool findCut(const vector<string>&, string);
  
  void updateMet(int syst=0);
  //  void treatMuons_Met(string syst="orig");
  double getPileupWeight(float);
  unordered_map<CUTS, vector<int>*, EnumHash> getArray();

  double getCRVal(string);
  void setupCR(string, double);

  ///// values /////

  TChain* BOOM;
  TTree* BAAM;
  TFile* infoFile;
  string filespace = "";
  double hPU[100];
  int version=0;

  Generated* _Gen;
  Electron* _Electron;
  Muon* _Muon;
  Taus* _Tau;
  Jet* _Jet;
  FatJet* _FatJet;
  Met* _MET;
  Histogramer histo;
  Histogramer syst_histo;
  Systematics systematics;
  JetScaleResolution jetScaleRes;
  PartStats genStat;

  unordered_map<string, PartStats> distats;
  unordered_map<string, FillVals*> fillInfo;
  unordered_map<string, double> genMap;
  unordered_map<CUTS, vector<int>*, EnumHash>* active_part;
  unordered_map<CUTS, vector<int>*, EnumHash> goodParts;
  vector<unordered_map<CUTS, vector<int>*, EnumHash>> syst_parts;
  unordered_map<CUTS, bool, EnumHash> need_cut;
  unordered_map<string,bool> gen_selection;

  vector<Particle*> allParticles;
  vector<string> syst_names;

  DepGraph neededCuts;

  static const unordered_map<string, CUTS> cut_num;
  static const unordered_map<CUTS, vector<CUTS>, EnumHash> adjList;

  vector<int>* trigPlace[nTrigReq];
  bool setTrigger = false;
  vector<string>* trigName[nTrigReq];
  vector<int> cuts_per, cuts_cumul;

  double maxIso, minIso;
  int leadIndex, maxCut, crbins=1;
  bool isData, CalculatePUSystematics, doSystematics;

  vector<int>* Trigger_decision = 0;
  vector<int>* Trigger_decisionV1 = 0;
  vector<string>* Trigger_names = 0;
  float nTruePU = 0;
  int bestVertices = 0;
  double gen_weight = 0;

  BTagCalibration calib = BTagCalibration("csvv1", "Pileup/btagging.csv");
  BTagCalibrationReader reader = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central");

  double rho =20.;

  const static vector<CUTS> genCuts;
  const static vector<CUTS> jetCuts;
  const static vector<CUTS> nonParticleCuts;
  double pu_weight, wgt, backup_wgt;
  unordered_map<int, GenFill*> genMaper;

  vector<CRTester*> testVec;
  int SignalRegion = -1;
  bool blinded = true;
  clock_t start_time;
  std::chrono::time_point<std::chrono::system_clock> start;

};



#endif

