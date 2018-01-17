#ifndef Particle_h
#define Particle_h

// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <regex>

#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>

#include "tokenizer.hpp"
#include "Cut_enum.h"

using namespace std;
typedef unsigned int uint;

struct PartStats {
  unordered_map<string,double> dmap;
  unordered_map<string,string> smap;
  unordered_map<string,pair<double,double> > pmap;
  //  unordered_map<string,bool> bmap;
  vector<string> bset;

  bool bfind(string cut) const {
    return find(bset.begin(), bset.end(), cut) != bset.end();
  }
  
};


enum class PType { Electron, Muon, Tau, Jet, FatJet, None};

class Particle {

public:
  Particle();
  Particle(TTree*, string, string, vector<string>);
  virtual ~Particle() {}

  virtual vector<CUTS> findExtraCuts() {return vector<CUTS>();}
  void init();
  void unBranch();
  double pt(uint) const;
  double eta(uint) const;
  double phi(uint) const;
  double energy(uint) const;
  virtual double charge(uint) const;
  TLorentzVector p4(uint) const;
  TLorentzVector& p4(uint);
  TLorentzVector RecoP4(uint) const;
  TLorentzVector& RecoP4(uint);

  uint size() const;
  vector<TLorentzVector>::iterator begin();
  vector<TLorentzVector>::iterator end();
  vector<TLorentzVector>::const_iterator begin() const;
  vector<TLorentzVector>::const_iterator end() const;

  bool needSyst(int) const;

  void addPtEtaPhiESyst(double, double, double, double, int);
  void addP4Syst(TLorentzVector, int);
  void setOrigReco();
  void setCurrentP(int);
  string getName() {return GenName;};

  bool findCut(const vector<string>&, string);
  
  PType type;
  unordered_map<string, PartStats> pstats;


protected:
  void getPartStats(string);
  TTree* BOOM;
  string GenName;
  const map<PType,CUTS> cutMap = {{PType::Electron, CUTS::eGElec}, {PType::Muon, CUTS::eGMuon},
				  {PType::Tau, CUTS::eGTau}};
  unordered_map<CUTS, string, EnumHash> jetNameMap = {
    {CUTS::eRJet1, "Jet1"},               {CUTS::eRJet2, "Jet2"},
    {CUTS::eRCenJet, "CentralJet"},      {CUTS::eRBJet, "BJet"},
    {CUTS::eR1stJet, "FirstLeadingJet"},  {CUTS::eR2ndJet, "SecondLeadingJet"},
    {CUTS::eRWjet, "WJet"}
  };

 private:
  vector<double>* mpt = 0;
  vector<double>* meta = 0;
  vector<double>* mphi = 0;
  vector<double>* menergy = 0;

  vector<TLorentzVector> Reco;
  vector<TLorentzVector> *cur_P;
  vector<string> syst_names;
  vector<vector<TLorentzVector>* > systVec;

  string activeSystematic;
};

class Photon : public Particle {
public:
  Photon();
  Photon(TTree*, string, vector<string>);

  vector<double>* et = 0;
  vector<double>* hoverE = 0;
  vector<double>* phoR = 0;
  vector<double>* sigmaIEtaIEta = 0;
  vector<double>* sigmaIPhiIPhi = 0;
  vector<double>* pfChIso = 0;
  vector<double>* pfPhoIso = 0;
  vector<double>* pfNeuIso = 0;
  vector<bool>*   eleVeto = 0;
  vector<bool>*   hasPixelSeed = 0;
};


/////////////////////////////////////////////////////////////////
class Generated : public Particle {

public:
  Generated();
  Generated(TTree*, string, vector<string>);

  vector<double>  *pdg_id = 0;
  vector<double>  *motherpdg_id = 0;
  vector<double>  *status = 0;
  vector<int>  *BmotherIndex = 0;

};

/////////////////////////////////////////////////////////////////////////
class Jet : public Particle {

public:
  Jet(TTree*, string, vector<string>);

  vector<CUTS> findExtraCuts();
  vector<CUTS> overlapCuts(CUTS);
  bool passedLooseJetID(int);
  
  vector<double>* neutralHadEnergyFraction = 0;
  vector<double>* neutralEmEmEnergyFraction = 0;
  vector<int>*    numberOfConstituents = 0;
  vector<double>* muonEnergyFraction = 0;
  vector<double>* chargedHadronEnergyFraction = 0;
  vector<int>*    chargedMultiplicity = 0;
  vector<double>* chargedEmEnergyFraction = 0;
  vector<int>*    partonFlavour = 0;
  vector<double>* bDiscriminator = 0;
  vector<double>* tau1 = 0;
  vector<double>* tau2 = 0;
  vector<double>* tau3 = 0;
  vector<double>* PrunedMass = 0;
  vector<double>* SoftDropMass = 0;

 protected:

};

class FatJet : public Particle {

public:
  FatJet(TTree*, string, vector<string>);

  vector<CUTS> findExtraCuts();
  vector<CUTS> overlapCuts(CUTS);

  vector<double>* tau1 = 0;
  vector<double>* tau2 = 0;
  vector<double>* tau3 = 0;
  vector<double>* PrunedMass = 0;
  vector<double>* SoftDropMass = 0;

};

class Lepton : public Particle {

public:
  Lepton(TTree*, string, string, vector<string>);

  vector<CUTS> findExtraCuts();

  double charge(uint)const;
  vector<double>* _charge = 0;


  virtual bool get_Iso(int, double, double) const {return false;}
};

class Electron : public Lepton {

public:
  Electron(TTree*, string, vector<string>);

  bool get_Iso(int, double, double) const;

  vector<int>     *isPassVeto = 0;
  vector<int>     *isPassLoose = 0;
  vector<int>     *isPassMedium = 0;
  vector<int>     *isPassTight = 0;
  vector<int>     *isPassHEEPId = 0;
  vector<double>  *isoChargedHadrons = 0;
  vector<double>  *isoNeutralHadrons = 0;
  vector<double>  *isoPhotons = 0;
  vector<double>  *isoPU = 0;
};



class Muon : public Lepton {

public:
  Muon(TTree*, string, vector<string>);

  bool get_Iso(int, double, double) const;

  vector<bool>* tight = 0;
  vector<bool>* soft = 0;
  vector<double>* isoCharged = 0;
  vector<double>* isoNeutralHadron = 0;
  vector<double>* isoPhoton = 0;
  vector<double>* isoPU = 0;
};

class Taus : public Lepton {

public:
  Taus(TTree*, string, vector<string>);

  //  void findExtraCuts();
  vector<CUTS> findExtraCuts();
  bool get_Iso(int, double, double) const;
  bool pass_against_Elec(CUTS, int);
  bool pass_against_Muon(CUTS, int);

  vector<int>     *decayModeFindingNewDMs = 0;
  vector<int>     *decayModeFinding = 0;
  vector<double>  *nProngs = 0;
  vector<int>  *decayMode = 0;
  pair<vector<int>*,vector<int>* > againstElectron = make_pair(nullptr,nullptr);
  pair<vector<int>*,vector<int>* > againstMuon = make_pair(nullptr,nullptr);
  pair<vector<int>*,vector<int>* > minIso = make_pair(nullptr,nullptr);
  pair<vector<int>*,vector<int>* > maxIso = make_pair(nullptr,nullptr);
  vector<double>  *leadChargedCandPt = 0;
  vector<double>  *leadChargedCandPtError = 0;
  vector<double>  *leadChargedCandValidHits = 0;
};



#endif
