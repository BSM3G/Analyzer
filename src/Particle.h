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

//using namespace std;
typedef unsigned int uint;

struct PartStats {
  std::unordered_map<std::string,double> dmap;
  std::unordered_map<std::string,std::string> smap;
  std::unordered_map<std::string,std::pair<double,double> > pmap;
  //  std::unordered_map<std::string,bool> bmap;
  std::vector<std::string> bset;

  bool bfind(std::string cut) const {
    return find(bset.begin(), bset.end(), cut) != bset.end();
  }
  
};


enum class PType { Electron, Muon, Tau, Jet, FatJet, None};


class Particle {

public:
  Particle();
  Particle(TTree*, std::string, std::string, std::vector<std::string>);
  virtual ~Particle() {}

  virtual std::vector<CUTS> findExtraCuts() {return std::vector<CUTS>();}
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
  std::vector<TLorentzVector>::iterator begin();
  std::vector<TLorentzVector>::iterator end();
  std::vector<TLorentzVector>::const_iterator begin() const;
  std::vector<TLorentzVector>::const_iterator end() const;

  bool needSyst(int) const;

  void addPtEtaPhiESyst(double, double, double, double, int);
  void addP4Syst(TLorentzVector, int);
  void setOrigReco();
  void setCurrentP(int);
  std::string getName() {return GenName;};

  bool findCut(const std::vector<std::string>&, std::string);
  
  PType type;
  std::unordered_map<std::string, PartStats> pstats;
  const std::map<PType,CUTS> cutMap = {{PType::Electron, CUTS::eGElec}, {PType::Muon, CUTS::eGMuon},
				  {PType::Tau, CUTS::eGTau}};


protected:
  void getPartStats(std::string);
  TTree* BOOM;
  std::string GenName;
  std::unordered_map<CUTS, std::string, EnumHash> jetNameMap = {
    {CUTS::eRJet1, "Jet1"},               {CUTS::eRJet2, "Jet2"},
    {CUTS::eRCenJet, "CentralJet"},      {CUTS::eRBJet, "BJet"},
    {CUTS::eR1stJet, "FirstLeadingJet"},  {CUTS::eR2ndJet, "SecondLeadingJet"},
    {CUTS::eRWjet, "WJet"}
  };

 private:
  std::vector<double>* mpt = 0;
  std::vector<double>* meta = 0;
  std::vector<double>* mphi = 0;
  std::vector<double>* menergy = 0;

  std::vector<TLorentzVector> Reco;
  std::vector<TLorentzVector> *cur_P;
  std::vector<std::string> syst_names;
  std::vector<std::vector<TLorentzVector>* > systVec;

  std::string activeSystematic;
};

class Photon : public Particle {
public:
  Photon();
  Photon(TTree*, std::string, std::vector<std::string>);

  std::vector<double>* et = 0;
  std::vector<double>* hoverE = 0;
  std::vector<double>* phoR = 0;
  std::vector<double>* sigmaIEtaIEta = 0;
  std::vector<double>* sigmaIPhiIPhi = 0;
  std::vector<double>* pfChIso = 0;
  std::vector<double>* pfPhoIso = 0;
  std::vector<double>* pfNeuIso = 0;
  std::vector<bool>*   eleVeto = 0;
  std::vector<bool>*   hasPixelSeed = 0;
};


/////////////////////////////////////////////////////////////////
class Generated : public Particle {

public:
  Generated();
  Generated(TTree*, std::string, std::vector<std::string>);

  std::vector<double>  *pdg_id = 0;
  std::vector<double>  *motherpdg_id = 0;
  std::vector<double>  *status = 0;
  std::vector<int>  *BmotherIndex = 0;
  std::vector<double> *numDaught = 0;     //04.22.18

};

/////////////////////////////////////////////////////////////////////////
class Jet : public Particle {

public:
  Jet(TTree*, std::string, std::vector<std::string>);

  std::vector<CUTS> findExtraCuts();
  std::vector<CUTS> overlapCuts(CUTS);
  bool passedLooseJetID(int);
  
  std::vector<double>* neutralHadEnergyFraction = 0;
  std::vector<double>* neutralEmEmEnergyFraction = 0;
  std::vector<int>*    numberOfConstituents = 0;
  std::vector<double>* muonEnergyFraction = 0;
  std::vector<double>* chargedHadronEnergyFraction = 0;
  std::vector<int>*    chargedMultiplicity = 0;
  std::vector<double>* chargedEmEnergyFraction = 0;
  std::vector<int>*    partonFlavour = 0;
  std::vector<double>* bDiscriminator = 0;
  std::vector<double>* tau1 = 0;
  std::vector<double>* tau2 = 0;
  std::vector<double>* tau3 = 0;
  std::vector<double>* PrunedMass = 0;
  std::vector<double>* SoftDropMass = 0;

 protected:

};

class FatJet : public Particle {

public:
  FatJet(TTree*, std::string, std::vector<std::string>);

  std::vector<CUTS> findExtraCuts();
  std::vector<CUTS> overlapCuts(CUTS);

  std::vector<double>* tau1 = 0;
  std::vector<double>* tau2 = 0;
  std::vector<double>* tau3 = 0;
  std::vector<double>* PrunedMass = 0;
  std::vector<double>* SoftDropMass = 0;

};

class Lepton : public Particle {

public:
  Lepton(TTree*, std::string, std::string, std::vector<std::string>);

  std::vector<CUTS> findExtraCuts();

  double charge(uint)const;
  std::vector<double>* _charge = 0;


  virtual bool get_Iso(int, double, double) const {return false;}
};

class Electron : public Lepton {

public:
  Electron(TTree*, std::string, std::vector<std::string>);

  bool get_Iso(int, double, double) const;

  std::vector<int>     *isPassVeto = 0;
  std::vector<int>     *isPassLoose = 0;
  std::vector<int>     *isPassMedium = 0;
  std::vector<int>     *isPassTight = 0;
  std::vector<int>     *isPassHEEPId = 0;
  std::vector<double>  *isoChargedHadrons = 0;
  std::vector<double>  *isoNeutralHadrons = 0;
  std::vector<double>  *isoPhotons = 0;
  std::vector<double>  *isoPU = 0;
};



class Muon : public Lepton {

public:
  Muon(TTree*, std::string, std::vector<std::string>);

  bool get_Iso(int, double, double) const;

  std::vector<bool>* tight = 0;
  std::vector<bool>* soft = 0;
  std::vector<double>* isoCharged = 0;
  std::vector<double>* isoNeutralHadron = 0;
  std::vector<double>* isoPhoton = 0;
  std::vector<double>* isoPU = 0;
};

class Taus : public Lepton {

public:
  Taus(TTree*, std::string, std::vector<std::string>);

  //  void findExtraCuts();
  std::vector<CUTS> findExtraCuts();
  bool get_Iso(int, double, double) const;
  bool pass_against_Elec(CUTS, int);
  bool pass_against_Muon(CUTS, int);

  std::vector<int>     *decayModeFindingNewDMs = 0;
  std::vector<int>     *decayModeFinding = 0;
  std::vector<double>  *nProngs = 0;
  std::vector<int>  *decayMode = 0;
  std::pair<std::vector<int>*,std::vector<int>* > againstElectron = std::make_pair(nullptr,nullptr);
  std::pair<std::vector<int>*,std::vector<int>* > againstMuon = std::make_pair(nullptr,nullptr);
  std::pair<std::vector<int>*,std::vector<int>* > minIso = std::make_pair(nullptr,nullptr);
  std::pair<std::vector<int>*,std::vector<int>* > maxIso = std::make_pair(nullptr,nullptr);
  std::vector<double>  *leadChargedCandPt = 0;
  std::vector<double>  *leadChargedCandPtError = 0;
  std::vector<double>  *leadChargedCandValidHits = 0;
  std::vector<double>  *leadChargedCandDz_pv = 0;
};



#endif
