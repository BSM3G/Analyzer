#include "SpechialAnalysis.h"
#include "HistClass.h"
#include <TFile.h>
#include <Compression.h>
#define BIG_NUM 46340

SpechialAnalysis::SpechialAnalysis(Analyzer* _a) {
  a=_a;
}

void SpechialAnalysis::init() {
  

  //this is just to store the histograms if they are to much to hold in memory
  // string safeFileName = "SpecialHistos.root";
  // file1               = new TFile(safeFileName.c_str(), "RECREATE");

  // number of events, saved in a histogram
  HistClass::CreateHistoUnchangedName("h_counters", 10, 0, 11, "N_{events}");

  for (unsigned int i = 0; i < 4; i++) {
    //str(boost::format("N_{%s}")%particleLatex[i] )
    HistClass::CreateHisto("num", particles[i].c_str(), 40, 0, 39, TString::Format("N_{%s}", particleSymbols[i].c_str()));
    HistClass::CreateHisto(3, "pt", particles[i].c_str(), 5000, 0, 5000, TString::Format("p_{T}^{%s} (GeV)", particleSymbols[i].c_str()));
    HistClass::CreateHisto(3, "eta", particles[i].c_str(), 80, -4, 4, TString::Format("#eta_{%s}", particleSymbols[i].c_str()));
    HistClass::CreateHisto(3, "phi", particles[i].c_str(), 40, -3.2, 3.2, TString::Format("#phi_{%s} (rad)", particleSymbols[i].c_str()));

    if (not a->isData) {
      HistClass::CreateHisto(1, "num_Gen", particles[i].c_str(), 40, 0, 39, TString::Format("N_{%s}", particleSymbols[i].c_str()));
      HistClass::CreateHisto(1, "pt_Gen", particles[i].c_str(), 5000, 0, 5000, TString::Format("p_{T}^{%s} (GeV)", particleSymbols[i].c_str()));
      HistClass::CreateHisto(1, "eta_Gen", particles[i].c_str(), 80, -4, 4, TString::Format("#eta_{%s}", particleSymbols[i].c_str()));
      HistClass::CreateHisto(1, "phi_Gen", particles[i].c_str(), 40, -3.2, 3.2, TString::Format("#phi_{%s} (rad)", particleSymbols[i].c_str()));
    }
  }
}

void SpechialAnalysis::begin_run() {
}

void SpechialAnalysis::analyze() {
  //just because it is convinient read in  the cuts:
  // const unordered_map<string,pair<int,int> >* cut_info = a->histo.get_cuts();
  // const vector<string>* cut_order = a->histo.get_cutorder();

  if(a->active_part->at(a->cut_num.at("NRecoTriggers1"))->size()==0 )
    return;

  int p1 = -1;
  int p2 = -1;
  if (a->active_part->at(CUTS::eDiTau)->size() == 1) {
    p1 = a->active_part->at(CUTS::eDiTau)->at(0) / BIG_NUM;
    p2 = a->active_part->at(CUTS::eDiTau)->at(0) % BIG_NUM;
  } else {
    return;
  }
  int j1      = -1;
  int j2      = -1;
  double mass = 0;
  for (auto it : *a->active_part->at(CUTS::eDiJet)) {
    int j1tmp = (it) / a->_Jet->size();
    int j2tmp = (it) % a->_Jet->size();
    if (a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "") > mass) {
      j1   = j1tmp;
      j2   = j2tmp;
      mass = a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "");
    }
  }
  if (p1 < 0 or p2 < 0 or j1 < 0 or j2 < 0)
    return;
    
  HistClass::Fill("Tau_num",a->active_part->at(CUTS::eRTau1)->size(),a->wgt);


}

void SpechialAnalysis::end_run() {
  
  TFile* outfile = new TFile(a->histo.outname.c_str(), "UPDATE", a->histo.outname.c_str(), ROOT::CompressionSettings(ROOT::kLZMA, 9));
  
  outfile->cd();
  outfile->mkdir("Spechial");
  outfile->cd("Spechial/");
  HistClass::WriteAll("Tau_");

  
  
}
