#include "Particle.h"

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

Particle::Particle(TTree* _BOOM, string _GenName, string filename) : BOOM(_BOOM), GenName(_GenName) {
  type = PType::None;
  getPartStats(filename);
  //set the pt to an empty vector if the branch does not exist
  if( _BOOM->GetListOfBranches()->FindObject((GenName+"_pt").c_str()) ==0){
    pt=new vector<double>();
  }else{
    SetBranch((GenName+"_pt").c_str(), pt);
    SetBranch((GenName+"_eta").c_str(), eta);
    SetBranch((GenName+"_phi").c_str(), phi);
    SetBranch((GenName+"_energy").c_str(), energy);
  }

}

void Particle::unBranch() {
  BOOM->SetBranchStatus((GenName+"*").c_str(), 0);
}


Photon::Photon(TTree* _BOOM, string filename) : Particle(_BOOM, "Photon", filename) {
  SetBranch("Photon_et", et);
  SetBranch("Photon_HoverE", hoverE);
  SetBranch("Photon_phoR9", phoR);
  SetBranch("Photon_SigmaIEtaIEta", sigmaIEtaIEta);
  SetBranch("Photon_SigmaIPhiIPhi", sigmaIPhiIPhi);
  SetBranch("Photon_PFChIso", pfChIso);
  SetBranch("Photon_PFPhoIso", pfPhoIso);
  SetBranch("Photon_PFNeuIso", pfNeuIso);
  SetBranch("Photon_EleVeto", eleVeto);
  SetBranch("Photon_hasPixelSeed", hasPixelSeed);
}

Generated::Generated(TTree* _BOOM, string filename) : Particle(_BOOM, "Gen", filename) {

  SetBranch("Gen_pdg_id", pdg_id);
  SetBranch("Gen_motherpdg_id", motherpdg_id);
  SetBranch("Gen_status", status);
  SetBranch("Gen_BmotherIndex", BmotherIndex);
}

Jet::Jet(TTree* _BOOM, string filename) : Particle(_BOOM, "Jet", filename) {
  type = PType::Jet;
  SetBranch("Jet_neutralHadEnergyFraction", neutralHadEnergyFraction);
  SetBranch("Jet_neutralEmEmEnergyFraction", neutralEmEmEnergyFraction);
  SetBranch("Jet_numberOfConstituents", numberOfConstituents);
  SetBranch("Jet_muonEnergyFraction", muonEnergyFraction);
  SetBranch("Jet_chargedHadronEnergyFraction", chargedHadronEnergyFraction);
  SetBranch("Jet_chargedMultiplicity", chargedMultiplicity);
  SetBranch("Jet_chargedEmEnergyFraction", chargedEmEnergyFraction);
  SetBranch("Jet_partonFlavour", partonFlavour);
  SetBranch("Jet_bDiscriminator_pfCISVV2", bDiscriminator);
}

void Jet::findExtraCuts() {
  if(pstats["Smear"].bmap.at("SmearTheJet")) {
    extraCuts.push_back(CUTS::eGMuon);
    extraCuts.push_back(CUTS::eGElec);
    extraCuts.push_back(CUTS::eGTau);
  }

}

vector<CUTS> Jet::overlapCuts(CUTS ePos) {
  string pName = jetNameMap.at(ePos);
  vector<CUTS> returnCuts;
  if(pstats.at(pName).bmap.at("RemoveOverlapWithMuon1s")) returnCuts.push_back(CUTS::eRMuon1);
  if(pstats[pName].bmap.at("RemoveOverlapWithMuon2s")) returnCuts.push_back(CUTS::eRMuon2);
  if(pstats[pName].bmap.at("RemoveOverlapWithElectron1s")) returnCuts.push_back(CUTS::eRElec1);
  if(pstats[pName].bmap.at("RemoveOverlapWithElectron2s")) returnCuts.push_back(CUTS::eRElec2);
  if(pstats[pName].bmap.at("RemoveOverlapWithTau1s")) returnCuts.push_back(CUTS::eRTau1);
  if(pstats[pName].bmap.at("RemoveOverlapWithTau2s")) returnCuts.push_back(CUTS::eRTau2);

  return returnCuts;
}


FatJet::FatJet(TTree* _BOOM, string filename) : Particle(_BOOM, "Jet_toptag", filename) {
  type = PType::FatJet;
  SetBranch("Jet_toptag_tau1", tau1);
  SetBranch("Jet_toptag_tau2", tau2);
  SetBranch("Jet_toptag_tau3", tau3);
  SetBranch("Jet_toptag_PrunedMass", PrunedMass);
  SetBranch("Jet_toptag_SoftDropMass", SoftDropMass);
}

void FatJet::findExtraCuts() {
  if(pstats["Smear"].bmap.at("SmearTheJet")) {
    extraCuts.push_back(CUTS::eGMuon);
    extraCuts.push_back(CUTS::eGElec);
    extraCuts.push_back(CUTS::eGTau);
  }
}

vector<CUTS> FatJet::overlapCuts(CUTS ePos) {
  string pName = jetNameMap.at(ePos);
  vector<CUTS> returnCuts;
  if(pstats.at(pName).bmap.at("RemoveOverlapWithMuon1s")) returnCuts.push_back(CUTS::eRMuon1);
  if(pstats[pName].bmap.at("RemoveOverlapWithMuon2s")) returnCuts.push_back(CUTS::eRMuon2);
  if(pstats[pName].bmap.at("RemoveOverlapWithElectron1s")) returnCuts.push_back(CUTS::eRElec1);
  if(pstats[pName].bmap.at("RemoveOverlapWithElectron2s")) returnCuts.push_back(CUTS::eRElec2);
  if(pstats[pName].bmap.at("RemoveOverlapWithTau1s")) returnCuts.push_back(CUTS::eRTau1);
  if(pstats[pName].bmap.at("RemoveOverlapWithTau2s")) returnCuts.push_back(CUTS::eRTau2);

  return returnCuts;
}


Lepton::Lepton(TTree* _BOOM, string GenName, string EndName) : Particle(_BOOM, GenName, EndName) {
  SetBranch((GenName+"_charge").c_str(), charge);
}

void Lepton::findExtraCuts() {
  if(pstats["Smear"].bmap.at("SmearTheParticle") || pstats["Smear"].bmap.at("MatchToGen")) {
    extraCuts.push_back(cutMap.at(type));
  }
}

Electron::Electron(TTree* _BOOM, string filename) : Lepton(_BOOM, "patElectron", filename) {
  type = PType::Electron;
  if(pstats["Elec1"].bmap["DoDiscrByIsolation"] || pstats["Elec2"].bmap["DoDiscrByIsolation"]) {
    SetBranch("patElectron_isoChargedHadrons", isoChargedHadrons);
    SetBranch("patElectron_isoNeutralHadrons", isoNeutralHadrons);
    SetBranch("patElectron_isoPhotons", isoPhotons);
    SetBranch("patElectron_isoPU", isoPU);
  }
  if(pstats["Elec1"].bmap["DoDiscrByVetoID"] || pstats["Elec2"].bmap["DoDiscrByVetoID"]) {
    SetBranch("patElectron_isPassVeto", isPassVeto);
  }
  if(pstats["Elec1"].bmap["DoDiscrByLooseID"] || pstats["Elec2"].bmap["DoDiscrByLooseID"]) {
    SetBranch("patElectron_isPassLoose", isPassLoose);
  }
  if(pstats["Elec1"].bmap["DoDiscrByMediumID"] || pstats["Elec2"].bmap["DoDiscrByMediumID"]) {
    SetBranch("patElectron_isPassMedium", isPassMedium);
  }
  if(pstats["Elec1"].bmap["DoDiscrByTightID"] || pstats["Elec2"].bmap["DoDiscrByTightID"]) {
    SetBranch("patElectron_isPassTight", isPassTight);
  }
  if(pstats["Elec1"].bmap["DoDiscrByHEEPID"] || pstats["Elec2"].bmap["DoDiscrByHEEPID"]) {
    SetBranch("patElectron_isPassHEEPId", isPassHEEPId);
  }
}

Muon::Muon(TTree* _BOOM, string filename) : Lepton(_BOOM, "Muon", filename) {
  type = PType::Muon;

  if(pstats["Muon1"].bmap["DoDiscrByTightID"] || pstats["Muon2"].bmap["DoDiscrByTightID"]) {
    SetBranch("Muon_tight", tight);
  }
  if(pstats["Muon1"].bmap["DoDiscrBySoftID"] || pstats["Muon2"].bmap["DoDiscrBySoftID"]) {
    SetBranch("Muon_soft", soft);
  }
  if(pstats["Muon1"].bmap["DoDiscrByIsolation"] || pstats["Muon2"].bmap["DoDiscrByIsolation"]) {
    SetBranch("Muon_isoCharged", isoCharged);
    SetBranch("Muon_isoNeutralHadron", isoNeutralHadron);
    SetBranch("Muon_isoPhoton", isoPhoton);
    SetBranch("Muon_isoPU", isoPU);
  }
}

///////fix against stuff
Taus::Taus(TTree* _BOOM, string filename) : Lepton(_BOOM, "Tau", filename) {
  type = PType::Tau;

  ////Electron discrimination
  if((pstats["Tau1"].bmap["DoDiscrAgainstElectron"] || pstats["Tau1"].bmap["SelectTausThatAreElectrons"]) &&
  (pstats["Tau2"].bmap["DoDiscrAgainstElectron"] || pstats["Tau2"].bmap["SelectTausThatAreElectrons"]) &&
  (pstats["Tau1"].smap["DiscrAgainstElectron"] == pstats["Tau2"].smap["DiscrAgainstElectron"]) ) {
    SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrAgainstElectron"]).c_str(), againstElectron.first);
    againstElectron.second = againstElectron.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrAgainstElectron"] || pstats["Tau1"].bmap["SelectTausThatAreElectrons"]) {
      SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrAgainstElectron"]).c_str(), againstElectron.first);
    }
    if(pstats["Tau2"].bmap["DoDiscrAgainstElectron"] || pstats["Tau2"].bmap["SelectTausThatAreElectrons"]) {
      SetBranch(("Tau_"+pstats["Tau2"].smap["DiscrAgainstElectron"]).c_str(), againstElectron.second);
    }
  }
  ////Muon discrimination
  if((pstats["Tau1"].bmap["DoDiscrAgainstMuon"] || pstats["Tau1"].bmap["SelectTausThatAreMuons"]) &&
  (pstats["Tau2"].bmap["DoDiscrAgainstMuon"] || pstats["Tau2"].bmap["SelectTausThatAreMuons"]) &&
  (pstats["Tau1"].smap["DiscrAgainstMuon"] == pstats["Tau2"].smap["DiscrAgainstMuon"]) ) {
    SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrAgainstMuon"]).c_str(), againstMuon.first);
    againstMuon.second = againstMuon.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrAgainstMuon"] || pstats["Tau1"].bmap["SelectTausThatAreMuons"]) {
      SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrAgainstMuon"]).c_str(), againstMuon.first);
    }
    if(pstats["Tau2"].bmap["DoDiscrAgainstMuon"] || pstats["Tau2"].bmap["SelectTausThatAreMuons"]) {
      SetBranch(("Tau_"+pstats["Tau2"].smap["DiscrAgainstMuon"]).c_str(), againstMuon.second);
    }
  }

  /////Isolation discrimination
  if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].bmap["DoDiscrByIsolation"] &&
  pstats["Tau1"].smap["DiscrByMaxIsolation"] == pstats["Tau2"].smap["DiscrByMaxIsolation"]) {
    SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrByMaxIsolation"]).c_str(), (maxIso.first));
    maxIso.second = maxIso.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrByIsolation"]) {
      SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrByMaxIsolation"]).c_str(), (maxIso.first));
    }
    if(pstats["Tau2"].bmap["DoDiscrByIsolation"]) {
      SetBranch(("Tau_"+pstats["Tau2"].smap["DiscrByMaxIsolation"]).c_str(), (maxIso.second));
    }
  }      ////min stuff
  if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].bmap["DoDiscrByIsolation"] &&
  pstats["Tau1"].smap["DiscrByMinIsolation"] == pstats["Tau2"].smap["DiscrByMinIsolation"] &&
  pstats["Tau1"].smap["DiscrByMinIsolation"] != "ZERO") {
    SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrByMinIsolation"]).c_str(), (minIso.first));
    minIso.second = minIso.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau1"].smap["DiscrByMinIsolation"] != "ZERO") {

      SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrByMinIsolation"]).c_str(), (minIso.first));
    }
    if(pstats["Tau2"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].smap["DiscrByMinIsolation"] != "ZERO") {
      SetBranch(("Tau_"+pstats["Tau2"].smap["DiscrByMaxIsolation"]).c_str(), (minIso.second));
    }
  }


  SetBranch("Tau_decayModeFindingNewDMs", decayModeFindingNewDMs);
  SetBranch("Tau_nProngs", nProngs);
  SetBranch("Tau_leadChargedCandPt", leadChargedCandPt);

}

void Taus::findExtraCuts() {
  Lepton::findExtraCuts();

  if(pstats["Tau1"].bmap.at("RemoveOverlapWithMuon1s") ||pstats["Tau2"].bmap.at("RemoveOverlapWithMuon1s"))
  extraCuts.push_back(CUTS::eRMuon1);
  if(pstats["Tau1"].bmap.at("RemoveOverlapWithMuon2s") ||pstats["Tau2"].bmap.at("RemoveOverlapWithMuon2s"))
  extraCuts.push_back(CUTS::eRMuon2);
  if(pstats["Tau1"].bmap.at("RemoveOverlapWithElectron1s") ||pstats["Tau2"].bmap.at("RemoveOverlapWithElectron1s"))
  extraCuts.push_back(CUTS::eRElec1);
  if(pstats["Tau1"].bmap.at("RemoveOverlapWithElectron2s") ||pstats["Tau2"].bmap.at("RemoveOverlapWithElectron2s"))
  extraCuts.push_back(CUTS::eRElec2);
}

void Particle::getPartStats(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    return;
  }

  vector<string> stemp;
  string group,line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      cout << "error in " << filename << "; no groups specified for data" << endl;
      exit(1);
    } else if(stemp.size() == 2) {

      if(stemp[1] == "1" || stemp[1] == "true" ) pstats[group].bmap[stemp[0]] = true;
      else if(stemp[1] == "0"  || stemp[1] == "false" ) pstats[group].bmap[stemp[0]]=false;

      else if(stemp[1].find_first_not_of("0123456789+-.") == string::npos) pstats[group].dmap[stemp[0]]=stod(stemp[1]);
      else pstats[group].smap[stemp[0]] = stemp[1];
    } else  pstats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
  }
  info_file.close();


}




bool Electron::get_Iso(int index, double min, double max) const {
  double maxIsoval = std::max(0.0, isoNeutralHadrons->at(index) + isoPhotons->at(index) - 0.5 * isoPU->at(index));
  double isoSum = (isoChargedHadrons->at(index) + maxIsoval) / smearP.at(index).Pt();
  return (isoSum >= min && isoSum < max);
}

bool Muon::get_Iso(int index, double min, double max) const {
  double maxIsoval = std::max(0.0, isoNeutralHadron->at(index) + isoPhoton->at(index) - 0.5 * isoPU->at(index));
  double isoSum = (isoCharged->at(index) + maxIsoval) / smearP.at(index).Pt();
  return (isoSum >= min && isoSum < max);
}

bool Taus::get_Iso(int index, double onetwo, double max) const {
  double maxIsoval = (onetwo == 1) ? maxIso.first->at(index) : maxIso.second->at(index);
  vector<int>* minIsotmp = (onetwo == 1) ? minIso.first : minIso.second;
  double minIsoval = (minIsotmp != 0) ? minIsotmp->at(index) : true;
  return (maxIsoval > 0.5 && minIsoval > 0.5);
}
