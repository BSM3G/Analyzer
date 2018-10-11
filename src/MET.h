#ifndef Met_h
#define Met_h

// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include "Particle.h"
#include "mt2/mt2_bisect.hh"


#include "tokenizer.hpp"
#include "Cut_enum.h"

//using namespace std;
typedef unsigned int uint;


class Met {

public:
  Met(){};
  Met(TTree*, std::string, std::string, std::vector<std::string>){};
  Met(TTree*, std::string, std::vector<std::string>, double);
  virtual ~Met() {}

  virtual std::vector<CUTS> findExtraCuts(){return std::vector<CUTS>();}
  void init();
  void unBranch();
  double pt() const;
  double px() const;
  double py() const;
  double eta() const;
  double phi() const;
  double energy() const;
  double HT() const;
  double MHT() const;
  double MHTphi() const;
  double MT2(TLorentzVector&, TLorentzVector&);
  TLorentzVector p4() const;
  TLorentzVector& p4();

  void addPtEtaPhiESyst(double, double, double, double, int);
  void addP4Syst(TLorentzVector, int);
  void setMT2Mass(double);
  void setCurrentP(int);
  std::string getName() {return GenName;};
  void update(PartStats&, Jet&, int);

  TLorentzVector Reco;
  TLorentzVector *cur_P;

  std::vector<TLorentzVector* > systVec;
  std::vector<double> systdeltaMEx;
  std::vector<double> systdeltaMEy;
  std::vector<double> syst_HT;
  std::vector<double> syst_MHT;
  std::vector<double> syst_MHTphi;


  int activeSystematic;

protected:
  TTree* BOOM;
  std::string GenName;
  std::vector<std::string> syst_names;
  double MT2mass;
  
  float m_pt;
  float m_phi;
  //note this is only for pt and phi
  double MetUnclUp[2] = {0, 0};
  double MetUnclDown[2] = {0, 0};
  int Unclup=-1;
  int Uncldown=-1;
  mt2_bisect::mt2 mt2_event;
  

};



#endif
