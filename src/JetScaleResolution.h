#ifndef JETSCALERESOLUTION
#define JETSCALERESOLUTION

#include "TRandom3.h"
#include <unordered_map>
#include <TLorentzVector.h>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
//#include "Particle.h"

using namespace std;

template<typename T> T stringtotype(string s)
{
    T i;
    istringstream(s) >> i;
    return(i);
}

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

vector<string> string_split(const string& in, const vector<string> splits = {" "});

//BINNER
class Bin
{
	private:
		double min_;
		double max_;
	public:
		Bin(double min, double max) : min_(min), max_(max) {}
		Bin(double val) : min_(val), max_(val) {}

		double min() const {return min_;}
		double max() const {return max_;}

inline bool operator<(const Bin& B){
        if(B.min() == B.max() && (this->min() <= B.min() && this->max() > B.min()))
	{
		return(false);
	}
	else if(this->min() == this->max() && (B.min() <= this->min() && B.max() > this->min()))
	{
		return(false);
	}
	return this->min() < B.min();
};

};

inline bool operator<(const Bin& A, const Bin& B)
{
	if(B.min() == B.max() && (A.min() <= B.min() && A.max() > B.min()))
	{
		return(false);
	}
	else if(A.min() == A.max() && (B.min() <= A.min() && B.max() > A.min()))
	{
		return(false);
	}
	return A.min() < B.min();
};



class JetScaleResolution{
    public:
        JetScaleResolution();
        JetScaleResolution(const string& scalefilename, const string& parttype, const string& resolutionfile, const string& sfresolutionfile);
        void InitScale(const string& filename, const string& type);
        void InitResolution(const string& resolutionfile, const string& sffile);
        double GetRes(const TLorentzVector& jet,const TLorentzVector& genjet, double rho, double sigmares);
        double GetScale(const TLorentzVector& jet, bool isBjet, double sigmascale);

    private:
        TH1D* Heta = nullptr;
        vector<TH1D*> HptsP;
        vector<TH1D*> HptsM;
        vector<TH1D*> HptsPqcd;
        vector<TH1D*> HptsMqcd;
        vector<TH1D*> HptsPb;
        vector<TH1D*> HptsMb;
        TGraph* hlE =nullptr;
        TGraph* hlB =nullptr;
        TGraph* hbE =nullptr;
        TGraph* hbB =nullptr;

        map< Bin, map<Bin, vector<double> > > resinfo;
        map< Bin, vector<double> > ressf;
};

#endif /*JETSCALERESOLUTION*/
