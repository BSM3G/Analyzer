#ifndef Histo_h
#define Histo_h

// system include files
#include <memory>

// user include files
#include <Math/VectorUtil.h>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <boost/unordered_map.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <regex>
#include "DataBinner.h"
#include "tokenizer.hpp"


using namespace std;

class Histogramer {

public:
  Histogramer();
  Histogramer(int, string, string, string, bool, vector<string>&);
  Histogramer(const Histogramer&);
  Histogramer(Histogramer&&);
  Histogramer& operator=(const Histogramer&);
  Histogramer& operator=(Histogramer&&);
  ~Histogramer();

  const unordered_map<string,pair<int,int>>* get_cuts() const {return &cuts;}
  const vector<string>* get_cutorder() const {return &cut_order;}
  const vector<string>* get_groups() const {return &data_order;}
  const vector<string>* get_folders() const {return &folders;}
  int get_maxfolder() const {return (folderToCutNum.back()+1);}

  void addVal(double, string, int, string, double);
  void addVal(double, double, string, int, string, double);
  void fill_histogram();
  void setControlRegions();

private:
  TFile * outfile;
  string outname;
  int NFolders;
  int Npdf;
  bool isData, CR=false;

  unordered_map<string, pair<int,int>> cuts;
  vector<string> cut_order;

  vector<string> folders;
  vector<int> folderToCutNum;

  boost::unordered_map<string, DataBinner*> data;
  vector<string> data_order;

  void read_hist(string);
  void read_cuts(string filename, vector<string>&);
  void fillCRFolderNames(string, int, bool, const vector<string>&);

  string extractHistname(string, string) const;
};

#endif
