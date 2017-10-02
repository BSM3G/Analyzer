#include "MET.h"
#include <algorithm>

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, true);  BOOM->SetBranchAddress(name, &variable);

//particle is a objet that stores multiple versions of the particle candidates
Met::Met(TTree* _BOOM, string _GenName,  vector<string> _syst_names) : BOOM(_BOOM), GenName(_GenName), syst_names(_syst_names) {

  SetBranch((GenName+"_px").c_str(), mMet[0]);
  SetBranch((GenName+"_py").c_str(), mMet[1]);
  SetBranch((GenName+"_pz").c_str(), mMet[2]);

  systdeltaMEx.resize(syst_names.size());
  systdeltaMEy.resize(syst_names.size());
  syst_HT.resize(syst_names.size());
  syst_MHT.resize(syst_names.size());
  syst_MHTphi.resize(syst_names.size());
  

  for( auto name : syst_names) {
    if(name == "orig") 
      systVec.push_back(new TLorentzVector);
    else if(name.find("Met")!=string::npos){
      systVec.push_back(new TLorentzVector);
    }
    else
      systVec.push_back(nullptr);
  }

  if( std::find(syst_names.begin(), syst_names.end(), "Met_Uncl_Up") != syst_names.end() && _BOOM->GetListOfBranches()->FindObject((GenName+"_UnclEnshiftedPtUp").c_str()) !=0){
    SetBranch((GenName+"_UnclEnshiftedPtUp").c_str(), MetUnclUp[0]);
    SetBranch((GenName+"_UnclEnshiftedPhiUp").c_str(), MetUnclUp[1]);
    Unclup = std::find(syst_names.begin(), syst_names.end(), "Met_Uncl_Up") -syst_names.begin();
  }
  if( std::find(syst_names.begin(), syst_names.end(), "Met_Uncl_Down") != syst_names.end() && _BOOM->GetListOfBranches()->FindObject((GenName+"_UnclEnshiftedPtDown").c_str()) !=0){
    SetBranch((GenName+"_UnclEnshiftedPtDown").c_str(), MetUnclDown[0]);
    SetBranch((GenName+"_UnclEnshiftedPhiDown").c_str(), MetUnclDown[1]);
    Uncldown = std::find(syst_names.begin(), syst_names.end(), "Met_Uncl_Down")- syst_names.begin();
  }

  activeSystematic=0;
}


void Met::addPtEtaPhiESyst(double ipt,double ieta, double iphi, double ienergy, int syst){
  systVec.at(syst)->SetPtEtaPhiE(ipt,ieta,iphi,ienergy);
}


void Met::addP4Syst(TLorentzVector mp4, int syst){
  systVec.at(syst)->SetPtEtaPhiE(mp4.Pt(),mp4.Eta(),mp4.Phi(),mp4.E());
}


void Met::init(){
  //cleanup of the particles
  //keep this if there is any ever some need for a unchanged met
  Reco.SetPxPyPzE(mMet[0],mMet[1],mMet[2],sqrt(pow(mMet[0],2) + pow(mMet[1],2)));
  for(int i=0; i < (int) syst_names.size(); i++) {
    if(i == Unclup) systVec.at(i)->SetPtEtaPhiE(MetUnclUp[0],0,MetUnclUp[1],MetUnclUp[0]);
    else if(i == Uncldown) systVec.at(i)->SetPtEtaPhiE(MetUnclDown[0],0,MetUnclDown[1],MetUnclDown[0]);
    else if(systVec.at(i) != nullptr) addP4Syst(Reco, i);
    
    fill(systdeltaMEx.begin(), systdeltaMEx.end(), 0);
    fill(systdeltaMEy.begin(), systdeltaMEy.end(), 0);
  }
  cur_P=&Reco;

  // syst_HT[activeSystematic]=0.;
  // syst_MHT[activeSystematic]=0.;
  // syst_MHTphi[activeSystematic]=99.;
}


void Met::update(PartStats& stats, Jet& jet, int syst=0){
  ///Calculates met from values from each file plus smearing and treating muons as neutrinos
  if(systVec.at(syst) == nullptr) return;
  double sumpxForMht=0;
  double sumpyForMht=0;
  double sumptForHt=0;

  int i=0;
  for(auto jetVec: jet) {
    bool add = true;
    if( (jetVec.Pt() < stats.dmap.at("JetPtForMhtAndHt")) ||
	(abs(jetVec.Eta()) > stats.dmap.at("JetEtaForMhtAndHt")) ||
	( stats.bfind("ApplyJetLooseIDforMhtAndHt") &&
	  !jet.passedLooseJetID(i) ) ) add = false;
    if(add) {
      sumpxForMht -= jetVec.Px();
      sumpyForMht -= jetVec.Py();
      sumptForHt  += jetVec.Pt();
    }
    i++;
  }
  syst_HT.at(syst)=sumptForHt;
  syst_MHT.at(syst)=sumpyForMht;
  syst_MHTphi.at(syst)=atan2(sumpyForMht,sumpxForMht);

  systVec.at(syst)->SetPxPyPzE(systVec.at(syst)->Px()+systdeltaMEx[syst], 
			       systVec.at(syst)->Py()+systdeltaMEy[syst], 
			       systVec.at(syst)->Pz(), 
			       TMath::Sqrt(pow(systVec.at(syst)->Px()+systdeltaMEx[syst],2) + pow(systVec.at(syst)->Py()+systdeltaMEy[syst],2)));

}

double Met::pt()const         {return cur_P->Pt();}
double Met::px()const         {return cur_P->Px();}
double Met::py()const         {return cur_P->Py();}
double Met::eta()const        {return cur_P->Eta();}
double Met::phi()const        {return cur_P->Phi();}
double Met::energy()const     {return cur_P->E();}


TLorentzVector Met::p4()const {return (*cur_P);}
TLorentzVector& Met::p4() {return *cur_P;}

double Met::HT() const {return syst_HT.at(activeSystematic);};
double Met::MHT() const {return syst_MHT.at(activeSystematic);};
double Met::MHTphi() const {return syst_MHTphi.at(activeSystematic);};


void Met::setCurrentP(int syst){
  if( systVec.at(syst) == nullptr) {
    cur_P = systVec.at(0);
    activeSystematic = 0;
  } else {
    cur_P = systVec[syst];
    activeSystematic=syst;
  }
}


void Met::unBranch() {
  BOOM->SetBranchStatus((GenName+"*").c_str(), 0);
}

