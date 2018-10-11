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
#include <TTree.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <regex>
#include "DataBinner.h"
#include "tokenizer.hpp"


//using namespace std;

class Histogramer {

public:
  Histogramer();
  Histogramer(int, std::string, std::string, std::string, bool, std::vector<std::string>&, const std::vector<std::string> &syst_unvertainties={});
  Histogramer(const Histogramer&);
  Histogramer(Histogramer&&);
  Histogramer& operator=(const Histogramer&);
  Histogramer& operator=(Histogramer&&);
  ~Histogramer();

  const std::unordered_map<std::string,std::pair<int,int>>* get_cuts() const {return &cuts;}
  const std::vector<std::string>* get_cutorder() const {return &cut_order;}
  const std::vector<std::string>* get_groups() const {return &data_order;}
  const std::vector<std::string>* get_folders() const {return &folders;}
  int get_maxfolder() const {return (folderToCutNum.back()+1);}

  void addVal(double, std::string, int, std::string, double);
  void addVal(double, double, std::string, int, std::string, double);
  void addEffiency(std::string,double,bool,int);
  void fill_histogram(TFile* _outfile, std::string subfolder="");
  void setControlRegions();
  void createTree(std::unordered_map< std::string , float >*, std::string);
  void fillTree(std::string);
  std::string outname;


private:
  TFile * outfile;
  
  int NFolders;
  int Npdf;
  bool isData, fillSingle=false;

  std::unordered_map<std::string, std::pair<int,int>> cuts;
  std::vector<std::string> cut_order;

  std::vector<std::string> folders;
  std::vector<int> folderToCutNum;

  std::unordered_map<std::string, DataBinner*> data;
  std::vector<std::string> data_order;
  std::unordered_map<std::string, TTree * > trees;

  void read_hist(std::string);
  void read_cuts(std::string filename, std::vector<std::string>&);
  void read_syst(const std::vector<std::string>& syst_uncertainties);
  void fillCRFolderNames(std::string, int, bool, const std::vector<std::string>&);

  std::string extractHistname(std::string, std::string) const;
};

#endif
