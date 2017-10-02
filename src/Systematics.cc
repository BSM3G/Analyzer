#include "Systematics.h"


Systematics::Systematics(){}

Systematics::Systematics(std::unordered_map<string, PartStats> const &distats){

}
Systematics::~Systematics(){}

void Systematics::init(){

}


void Systematics::shiftParticle(Particle& jet, TLorentzVector recJet, double const& ratio, double& dPx, double& dPy, int syst){

   //add the shifted part up
   dPx+=recJet.Px()*(ratio-1);
   dPy+=recJet.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recJet*=ratio;
   jet.addP4Syst(recJet, syst);
   return;
}

void Systematics::shiftLepton(Lepton& lepton, TLorentzVector recoLep, TLorentzVector genLep, double& dPx, double& dPy, int syst){
  if (genLep == TLorentzVector(0,0,0,0)) {
    lepton.addP4Syst(recoLep, syst);
    return;
  }
  double ratio = (genLep.Pt()*scale) + (recoLep.Pt() - genLep.Pt())*resolution;
   //add the shifted part up
   dPx+=recoLep.Px()*(ratio-1);
   dPy+=recoLep.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recoLep*=ratio;
   lepton.addP4Syst(recoLep, syst);
   return;
}


void Systematics::loadScaleRes(const PartStats& smear, const PartStats& syst, string syst_name) {
  scale = 1;
  resolution = 1;
  if(smear.bfind("SmearTheParticle")) {
    scale = smear.dmap.at("PtScaleOffset");
    resolution = smear.dmap.at("PtScaleOffset");
  } 
  if(syst_name.find("_Res_")) {
    resolution = syst_name.find("_Up") ? 1 + syst.dmap.at("res") : 1 - syst.dmap.at("res");
  } else if(syst_name.find("_Res_")) {
    scale = syst_name.find("_Up") ? 1+syst.dmap.at("scale") : 1- syst.dmap.at("scale");
  }
}

