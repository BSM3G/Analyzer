#include "CRTest.h"
#define BIG_NUM 46340


bool CRTester::partPassBoth(Analyzer* analyzer) {
  int diffsize = 0;
  CUTS ePart1;
  CUTS ePart2;
  if(partName == "Muon1Muon2") {
    ePart1 = CUTS::eRMuon1;
    ePart2 = CUTS::eRMuon2;
    diffsize = analyzer->_Muon->smearP.size();      
  } else if(partName == "Electron1Electron2") {
    ePart1 = CUTS::eRElec1;
    ePart2 = CUTS::eRElec2;
    diffsize = analyzer->_Electron->smearP.size();
  } else if(partName == "Tau1Tau2") {
    ePart1 = CUTS::eRTau1;
    ePart2 = CUTS::eRTau2;
    diffsize = analyzer->_Tau->smearP.size();
  }
  vector<int>* part1 = analyzer->goodParts[ePart1];
  vector<int>* part2 = analyzer->goodParts[ePart2];
  vector<int> diff(diffsize);

  vector<int>::iterator it = set_symmetric_difference(part1->begin(), part1->end(), part2->begin(), part2->end(), diff.begin());
  diff.resize(it - diff.begin());

  return (diff.size() == 0);
}


CRTester::CRTester(FillVals* _info, string var, double val, string _name) : info(_info), variable(var), cutVal(val), partName(_name) {}

bool CRTester::test(Analyzer* analyzer) {

  bool pass = true;

  if(variable == "Pass") {
    return ((int)analyzer->goodParts[info->ePos]->size() >= cutVal);
  } else if(info->type == FILLER::Single) {
    if(variable == "PassBoth") return partPassBoth(analyzer);

    for(auto index: *analyzer->getList(info->ePos)) {
  
      TLorentzVector part = info->part->smearP.at(index);
      if(variable == "Eta") pass = pass && part.Eta() > cutVal;
      else if(variable == "Pt") pass = pass && part.Pt() > cutVal;
      else if(variable == "Energy") pass = pass && part.Energy() > cutVal;

    }	  
  } else if(info->type == FILLER::Dipart) {
    for(auto index: *analyzer->getList(info->ePos)) {

      TLorentzVector part1 = info->part->smearP.at(index / BIG_NUM);
      TLorentzVector part2 = info->part2->smearP.at(index % BIG_NUM);
      if(variable == "DeltaR") pass = pass && part1.DeltaR(part2) > cutVal;
      else if(variable == "DeltaPtDivSumPt") pass = pass && ((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()) > cutVal);
      else if(variable == "DeltaPt") pass = pass && ((part1.Pt() - part2.Pt()) > cutVal);
      else if(variable == "Zeta") pass = pass && (analyzer->getZeta(part1, part2, partName) > cutVal);
      else if(variable == "CosDphi") pass = pass && absnormPhi( part1.Phi() - part2.Phi()) > cutVal;
      else if(variable == "Mass") pass = pass && analyzer->getMass(part1, part2, partName) > cutVal;
      else if(variable == "DeltaEta") pass = pass && (abs(part1.Eta() - part2.Eta()) > cutVal);
      else if(variable == "DeltaPhi") pass = pass && (abs(part1.Phi() - part2.Phi()) > cutVal);
      else if(variable == "OSEta") pass = pass && (part1.Eta() * part2.Eta() > cutVal);
    }   
  } else {

    if(variable == "Met") pass = (analyzer->getMet() > cutVal);
    else if(variable == "Ht") pass = (analyzer->getHT() > cutVal);
    else if(variable == "Mht") pass = (analyzer->getMHT() > cutVal);
  }


  return pass;
}
