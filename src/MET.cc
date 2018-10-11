#include "MET.h"
#include <algorithm>

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, true);  BOOM->SetBranchAddress(name, &variable);

//particle is a objet that stores multiple versions of the particle candidates
Met::Met(TTree* _BOOM, std::string _GenName,  std::vector<std::string> _syst_names, double _MT2mass) : BOOM(_BOOM), GenName(_GenName), syst_names(_syst_names), MT2mass(_MT2mass)  {

  SetBranch((GenName+"_pt").c_str(), m_pt);
  SetBranch((GenName+"_phi").c_str(), m_phi);

  systdeltaMEx.resize(syst_names.size());
  systdeltaMEy.resize(syst_names.size());
  syst_HT.resize(syst_names.size());
  syst_MHT.resize(syst_names.size());
  syst_MHTphi.resize(syst_names.size());
  

  for( auto name : syst_names) {
    if(name == "orig") 
      systVec.push_back(new TLorentzVector);
    else if(name.find("Met")!=std::string::npos){
      systVec.push_back(new TLorentzVector);
    }
    else if(name.find("weight")!=std::string::npos){
      systVec.push_back(nullptr);
    }else if(name.find("Tau_qcd")!=std::string::npos){
      systVec.push_back(nullptr);
    }else
      systVec.push_back(new TLorentzVector);
  }
  
  if( std::find(syst_names.begin(), syst_names.end(), "MetUncl_Up") != syst_names.end() && _BOOM->GetListOfBranches()->FindObject((GenName+"_MetUnclustEnUpDeltaX").c_str()) !=0){
    SetBranch((GenName+"_MetUnclustEnUpDeltaX").c_str(), MetUnclUp[0]);
    SetBranch((GenName+"_MetUnclustEnUpDeltaY").c_str(), MetUnclUp[1]);
    Unclup = std::find(syst_names.begin(), syst_names.end(), "MetUncl_Up") -syst_names.begin();
  }
  if( std::find(syst_names.begin(), syst_names.end(), "MetUncl_Down") != syst_names.end() && _BOOM->GetListOfBranches()->FindObject((GenName+"_MetUnclustEnUpDeltaY").c_str()) !=0){
    SetBranch((GenName+"_MetUnclustEnUpDeltaX").c_str(), MetUnclDown[0]);
    SetBranch((GenName+"_MetUnclustEnUpDeltaY").c_str(), MetUnclDown[1]);
    Uncldown = std::find(syst_names.begin(), syst_names.end(), "MetUncl_Down")- syst_names.begin();
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
  Reco.SetPtEtaPhiM(m_pt,0,m_phi,m_pt);
  
  for(int i=0; i < (int) syst_names.size(); i++) {
    
    if(i == Unclup) systVec.at(i)->SetPxPyPzE(MetUnclUp[0]+Reco.Px(),MetUnclUp[1]+Reco.Py(),0,sqrt(pow(MetUnclUp[0]+Reco.Px(),2)+pow(MetUnclUp[1]+Reco.Py(),2)));
    else if(i == Uncldown) systVec.at(i)->SetPxPyPzE(MetUnclDown[0]+Reco.Px(),MetUnclDown[1]+Reco.Py(),0,sqrt(pow(MetUnclDown[0]+Reco.Px(),2)+pow(MetUnclDown[1]+Reco.Py(),2)));
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
  syst_MHT.at(syst)= sqrt( pow(sumpxForMht,2.0) + pow(sumpyForMht,2.0) );
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

double Met::MT2(TLorentzVector& pa, TLorentzVector& pb){
  mt2_event.set_momenta( pa, pb, px(), py());
  mt2_event.set_mn( MT2mass );
  return mt2_event.get_mt2();
}

void Met::setMT2Mass(double mass){
  MT2mass=mass;
}

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

