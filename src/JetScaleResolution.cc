#include "JetScaleResolution.h"
#include <sstream>
#include <cmath>



JetScaleResolution::JetScaleResolution(){
}

JetScaleResolution::JetScaleResolution(const string& scalefilename, const string& parttype, const string& resolutionfile, const string& sfresolutionfile){
    InitScale(scalefilename, parttype);
    InitResolution(resolutionfile, sfresolutionfile);
}

void JetScaleResolution::InitScale(const string& filename, const string& type)
{
    vector<TH1D*>* ErrP = &HptsP;
    vector<TH1D*>* ErrM = &HptsM;

    //TDirectory* dir = gDirectory;
    //gROOT->cd();
    fstream fs(filename.c_str(), fstream::in);
    if(!fs.good()){
        cout<<"Jet file "<<filename<<" does not exist!"<<endl;
        exit(2);
    }
    string line;
    int selected = -1;
    vector<double> etabins;
    while(!fs.eof())
    {
      getline(fs, line);
      if(line.size() == 0) continue;
      selected++;
      if(selected>0){
        //cout<<line<<endl;
        vector<string> vals = string_split(line);
        if(selected == 0)
        {
            etabins.push_back(stringtotype<double>(vals[0]));
        }
        etabins.push_back(stringtotype<double>(vals[1]));
        vector<double> pts;
        vector<double> errm;
        vector<double> errp;
        for(size_t p = 3; p < vals.size() ; p+=3)
        {
            if(vals[p].size() == 0) break;
            pts.push_back(stringtotype<double>(vals[p]));
            errm.push_back(stringtotype<double>(vals[p+1]));
            errp.push_back(stringtotype<double>(vals[p+2]));
            //cout << p << " " << pts.back() << " " << errm.back() << " " << errp.back() << endl;
        }
        stringstream name;
        name << "jeterror_" << type << "_" << selected-2;
        ErrM->push_back(new TH1D((name.str()+"M").c_str(), (name.str()+"M").c_str(), pts.size()-1, pts.data()));
        ErrP->push_back(new TH1D((name.str()+"P").c_str(), (name.str()+"P").c_str(), pts.size()-1, pts.data()));
        for(size_t b = 0 ; b < pts.size()-1 ; ++b)
        {
            ErrM->back()->SetBinContent(b+1, errm[b]);
            ErrP->back()->SetBinContent(b+1, errp[b]);
        }
      }
    }
    Heta = new TH1D(("etabins_"+type).c_str(), ("etabins_"+type).c_str(), etabins.size()-1, etabins.data());
    if(etabins.size() == 0)
    {
        cerr << "ERROR - Jetscaler.InitScale: Could not find " << type << " in file " << filename << endl;
    }
    fs.close();
    //dir->cd();
}


void JetScaleResolution::InitResolution(const string& resolutionfile, const string& sffile)
{
    fstream fs(resolutionfile.c_str(), fstream::in);
    string line;
    while(!fs.eof())
    {
        getline(fs, line);
        if(line.size() == 0 || line[0] == '{') continue;
        vector<string> vals = string_split(line, {" ", "\t"});
        //for(string s:vals) cout << "|" << s << "|";
        //cout << endl;
        //cout << vals.size() << " " << vals[0] << " " << vals[1] << endl;
        if(vals.size() != 11) continue;

        Bin eta(stringtotype<double>(vals[0]), stringtotype<double>(vals[1]));
        Bin rho(stringtotype<double>(vals[2]), stringtotype<double>(vals[3]));
        //cout << stringtotype<double>(vals[0]) << " " << stringtotype<double>(vals[1]) << " " <<stringtotype<double>(vals[2]) << " " <<stringtotype<double>(vals[3]) << endl;
        resinfo[eta][rho] = {stringtotype<double>(vals[7]), stringtotype<double>(vals[8]), stringtotype<double>(vals[9]), stringtotype<double>(vals[10])};
    }
    fs.close();
    //cout << "finished "<< resinfo.size() << endl;
    fstream fsc(sffile.c_str(), fstream::in);
    while(!fsc.eof())
    {
        getline(fsc, line);
        if(line.size() == 0 || line[0] == '{') continue;
        vector<string> vals = string_split(line, {" ", "\t"});
        if(vals.size() != 6) continue;

        Bin eta(stringtotype<double>(vals[0]), stringtotype<double>(vals[1]));
        ressf[eta] = {stringtotype<double>(vals[3]), stringtotype<double>(vals[4]), stringtotype<double>(vals[5])};
    }
    fsc.close();
}



double JetScaleResolution::GetRes(const TLorentzVector& jet,const TLorentzVector& genjet, double rho, double sigmares)
{
    double rescor = 1.;
    if(rho > 44) {rho = 44;}
    if(abs(jet.Eta()) >= 4.7) {return 1.;}
    //cout << jet.Eta() << " " << rho << endl;
    const vector<double>& par = resinfo[jet.Eta()][rho];
    double x = jet.Pt();
    double resolution = sqrt(par[0]*abs(par[0])/(x*x)+par[1]*par[1]*pow(x,par[3])+par[2]*par[2]);
    //cout << jet.Eta() << " " << rho << " - " << resolution << ": " << par[0] << " " << par[1]<< " " << par[2] << endl;

    const vector<double>& sfs = ressf[jet.Eta()];
    if(sfs.size()==0){
      return(1.);
    }
    double s = sfs[0];
    //cout << sfs[0] << " " << sfs[1] << " " << sfs[2] << endl;
    if(sigmares <= 0) {s = sfs[0] + sigmares*(sfs[0]-sfs[1]);}
    if(sigmares > 0) {s = sfs[0] + sigmares*(sfs[2]-sfs[0]);}


    if(genjet == TLorentzVector(0,0,0,0) )
    {
        rescor +=  (s-1)*(jet.Pt()-genjet.Pt())/jet.Pt();
    }
    else
    {
        rescor += gRandom->Gaus(0., resolution*sqrt(s*s-1.));
    }
    return(max({0., rescor}));
}

double JetScaleResolution::GetScale(const TLorentzVector& jet, bool isBjet, double sigmascale)
{
    double sf = 1.;
    if(abs(jet.Eta()) >= 5.4) {return 1.;}

//  //cout << jet.Eta() << " " << rho << endl;
//  //cout << (resinfo.find(jet.Eta()) - resinfo.begin()) << endl;
//  const vector<double>& par = resinfo[jet.Eta()][rho];
//  double x = jet.Pt();
//  double resolution = sqrt(par[0]*abs(par[0])/(x*x)+par[1]*par[1]*pow(x,par[3])+par[2]*par[2]);
//  //cout << jet.Eta() << " " << rho << " - " << resolution << ": " << par[0] << " " << par[1]<< " " << par[2]<< " " << par[3] << endl;
//  const vector<double>& sfs = ressf[jet.Eta()];
//  double s = sfs[0];
//  //cout << sfs[0] << " " << sfs[1] << " " << sfs[2] << endl;
//  if(sigmares <= 0) {s = sfs[0] + sigmares*(sfs[0]-sfs[1]);}
//  if(sigmares > 0) {s = sfs[0] + sigmares*(sfs[2]-sfs[0]);}
//  sf = gRandom->Gaus(1., resolution*sqrt(s*s-1.));

    //MC specific correction
    //if(hlB != nullptr)
    //{
        //double mccorr = 1.;
        //if(isBjet)
        //{
            //if(abs(jet.Eta()) < 1.5)
            //{
                //mccorr = hbB->Eval(jet.Pt());
            //}
            //else
            //{
                //mccorr = hbE->Eval(jet.Pt());
            //}
        //}
        //else
        //{
            //if(abs(jet.Eta()) < 1.5)
            //{
                //mccorr = hlB->Eval(jet.Pt());
            //}
            //else
            //{
                //mccorr = hlE->Eval(jet.Pt());
            //}
        //}
        //sf += mccorr;
    //}

    int etabin = Heta->FindFixBin(jet.Eta());
    int ptbin = HptsP[etabin]->FindFixBin(jet.Pt());
    if(sigmascale >= 0)
    {
        if(isBjet) return(sf + sigmascale*sqrt(pow(HptsP[etabin]->GetBinContent(ptbin), 2) + pow(HptsPb[etabin]->GetBinContent(ptbin), 2)));
        return(sf + sigmascale*HptsP[etabin]->GetBinContent(ptbin));
    }
    else
    {
        if(isBjet) return(sf + sigmascale*sqrt(pow(HptsM[etabin]->GetBinContent(ptbin), 2) + pow(HptsMb[etabin]->GetBinContent(ptbin), 2)));
        return(sf + sigmascale*(HptsM[etabin]->GetBinContent(ptbin)));
    }
}

vector<string> string_split(const string& in, const vector<string> splits)
{
    vector<pair<size_t, size_t> > positions;
    positions.push_back(pair<size_t, size_t>(0, 0));
    for(size_t s = 0 ; s < splits.size() ; ++s)
    {
        size_t lastpos = 0;
        while(lastpos < in.size())
        {
            lastpos = in.find(splits[s], lastpos);
            if(lastpos == string::npos)
            {
                break;
            }
            else
            {
                positions.push_back(pair<size_t, size_t>(lastpos, splits[s].size()));
                //lastpos += splits[s].size()+1;
                lastpos += splits[s].size();
            }
        }

    }
    positions.push_back(pair<size_t, size_t>(in.size(), 0));
    sort(positions.begin(), positions.end(), [](const pair<size_t, size_t>& A, const pair<size_t, size_t>& B){return A.first < B.first;});
    vector<string> result;
    for(size_t p = 0 ; p < positions.size()-1 ; ++p)
    {
        size_t begin = positions[p].first + positions[p].second;
        size_t end = positions[p+1].first;
        if(end != begin)result.push_back(in.substr(begin, end-begin));
    }
    return result;
}


