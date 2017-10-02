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


#include "tokenizer.hpp"
#include "Cut_enum.h"

using namespace std;
typedef unsigned int uint;


class Met {

public:
  Met();
  Met(TTree*, string, vector<string>);
  virtual ~Met() {}

  virtual void findExtraCuts() {}
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
  TLorentzVector p4() const;
  TLorentzVector& p4();

  void addPtEtaPhiESyst(double, double, double, double, int);
  void addP4Syst(TLorentzVector, int);
  void setCurrentP(int);
  string getName() {return GenName;};
  void update(PartStats&, Jet&, int);

  TLorentzVector Reco;
  TLorentzVector *cur_P;

  vector<TLorentzVector* > systVec;
  vector<double> systdeltaMEx;
  vector<double> systdeltaMEy;
  vector<double> syst_HT;
  vector<double> syst_MHT;
  vector<double> syst_MHTphi;


  int activeSystematic;

protected:
  TTree* BOOM;
  string GenName;
  double mMet[3] = {0, 0, 0};
  //note this is only for pt and phi
  double MetUnclUp[2] = {0, 0};
  double MetUnclDown[2] = {0, 0};
  vector<string> syst_names;
  int Unclup=-1;
  int Uncldown=-1;

};



#endif
