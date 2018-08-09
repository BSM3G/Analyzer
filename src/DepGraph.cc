#include "DepGraph.h"

//using namespace boost;
//using namespace std;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> mygraph;
#define cutint(x) static_cast<int>(x)
#define intcut(x) static_cast<CUTS>(x)

DepGraph::DepGraph() {
  
  add_edge(cutint(CUTS::eGTau), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGTop), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGElec), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGMuon), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGZ), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGW), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGHiggs), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGJet), cutint(CUTS::eGen), g);
  
  add_edge(cutint(CUTS::eMuon1Tau1), cutint(CUTS::eRMuon1), g);
  add_edge(cutint(CUTS::eMuon1Tau1), cutint(CUTS::eRTau1), g);
  add_edge(cutint(CUTS::eMuon1Tau2), cutint(CUTS::eRMuon1), g);
  add_edge(cutint(CUTS::eMuon1Tau2), cutint(CUTS::eRTau2), g);
  add_edge(cutint(CUTS::eMuon2Tau1), cutint(CUTS::eRMuon2), g);
  add_edge(cutint(CUTS::eMuon2Tau1), cutint(CUTS::eRTau1), g);
  add_edge(cutint(CUTS::eMuon2Tau2), cutint(CUTS::eRMuon2), g);
  add_edge(cutint(CUTS::eMuon2Tau2), cutint(CUTS::eRTau2), g);

  add_edge(cutint(CUTS::eElec1Tau1), cutint(CUTS::eRElec1), g);
  add_edge(cutint(CUTS::eElec1Tau1), cutint(CUTS::eRTau1), g);
  add_edge(cutint(CUTS::eElec1Tau2), cutint(CUTS::eRElec1), g);
  add_edge(cutint(CUTS::eElec1Tau2), cutint(CUTS::eRTau2), g);
  add_edge(cutint(CUTS::eElec2Tau1), cutint(CUTS::eRElec2), g);
  add_edge(cutint(CUTS::eElec2Tau1), cutint(CUTS::eRTau1), g);
  add_edge(cutint(CUTS::eElec2Tau2), cutint(CUTS::eRElec2), g);
  add_edge(cutint(CUTS::eElec2Tau2), cutint(CUTS::eRTau2), g);

  add_edge(cutint(CUTS::eMuon1Elec1), cutint(CUTS::eRMuon1), g);
  add_edge(cutint(CUTS::eMuon1Elec1), cutint(CUTS::eRElec1), g);
  add_edge(cutint(CUTS::eMuon1Elec2), cutint(CUTS::eRMuon1), g);
  add_edge(cutint(CUTS::eMuon1Elec2), cutint(CUTS::eRElec2), g);
  add_edge(cutint(CUTS::eMuon2Elec1), cutint(CUTS::eRMuon2), g);
  add_edge(cutint(CUTS::eMuon2Elec1), cutint(CUTS::eRElec1), g);
  add_edge(cutint(CUTS::eMuon2Elec2), cutint(CUTS::eRMuon2), g);
  add_edge(cutint(CUTS::eMuon2Elec2), cutint(CUTS::eRElec2), g);

  add_edge(cutint(CUTS::eDiElec), cutint(CUTS::eRElec1), g);
  add_edge(cutint(CUTS::eDiElec), cutint(CUTS::eRElec2), g);
  add_edge(cutint(CUTS::eDiElec), cutint(CUTS::eR1stJet), g);
  add_edge(cutint(CUTS::eDiElec), cutint(CUTS::eR2ndJet), g);

  add_edge(cutint(CUTS::eDiMuon), cutint(CUTS::eRMuon1), g);
  add_edge(cutint(CUTS::eDiMuon), cutint(CUTS::eRMuon2), g);
  add_edge(cutint(CUTS::eDiMuon), cutint(CUTS::eR1stJet), g);
  add_edge(cutint(CUTS::eDiMuon), cutint(CUTS::eR2ndJet), g);

  add_edge(cutint(CUTS::eDiTau), cutint(CUTS::eRTau1), g);
  add_edge(cutint(CUTS::eDiTau), cutint(CUTS::eRTau2), g);
  add_edge(cutint(CUTS::eDiTau), cutint(CUTS::eR1stJet), g);
  add_edge(cutint(CUTS::eDiTau), cutint(CUTS::eR2ndJet), g);

  add_edge(cutint(CUTS::eDiJet), cutint(CUTS::eRJet1), g);
  add_edge(cutint(CUTS::eDiJet), cutint(CUTS::eRJet2), g);
  
  add_edge(cutint(CUTS::eElec1Jet1), cutint(CUTS::eRElec1), g);
  add_edge(cutint(CUTS::eElec1Jet1), cutint(CUTS::eRJet1), g);
  add_edge(cutint(CUTS::eElec1Jet2), cutint(CUTS::eRElec1), g);
  add_edge(cutint(CUTS::eElec1Jet2), cutint(CUTS::eRJet2), g);
  add_edge(cutint(CUTS::eElec2Jet1), cutint(CUTS::eRElec2), g);
  add_edge(cutint(CUTS::eElec2Jet1), cutint(CUTS::eRJet1), g);
  add_edge(cutint(CUTS::eElec2Jet2), cutint(CUTS::eRElec2), g);
  add_edge(cutint(CUTS::eElec2Jet2), cutint(CUTS::eRJet2), g);
  
  add_edge(cutint(CUTS::eSusyCom), cutint(CUTS::eR1stJet), g);
  add_edge(cutint(CUTS::eSusyCom), cutint(CUTS::eR2ndJet), g);

}

    
void DepGraph::dfs(int vertex) {
  mygraph::adjacency_iterator it, end;
  neededCuts.insert(vertex);

  for(tie(it, end) = adjacent_vertices(vertex, g); it != end; ++it) {
    if(neededCuts.find(*it) != neededCuts.end()) continue;
    dfs(*it);
  }
  return;
}


void DepGraph::loadCuts(std::vector<CUTS> cutVec) {
  for(auto cut: cutVec) {
    int icut = cutint(cut);
    if(neededCuts.find(icut) != neededCuts.end()) continue;
    
    dfs(icut);
  }
}

void DepGraph::loadCuts(CUTS cut) {
  int icut = cutint(cut);
  if(neededCuts.find(icut) != neededCuts.end()) return;

  dfs(icut);
}

bool DepGraph::isPresent(CUTS cut) {
  return (neededCuts.find(cutint(cut)) != neededCuts.end());
}

std::unordered_set<int> DepGraph::getCuts() {
  return neededCuts;
}
