
/** \namespace HistClass
 *
 * \brief Functions to have easy histogram handling
 *
 * In this namespace different functions to interact with the histogram
 * map are included, to create, fill and write histograms in a convinient
 * way.
 */

#ifndef HistClass_h
#define HistClass_h

#include <stdarg.h>
#include <iostream>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TString.h"
#include "TNtupleD.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "boost/format.hpp"

/** To avoid compiler problems, we tell gcc to ignore any unused function error
 */

#ifdef __GNUC__
#define SUPPRESS_NOT_USED_WARN __attribute__ ((unused))
#else
#define SUPPRESS_NOT_USED_WARN
#endif

namespace HistClass {
    static std::unordered_map<std::string, TH1D * > histo; /*!< Map of a string and a TH1D histogram, for easy 1D histogram handling. */
    static std::unordered_map<std::string, TH2D * > histo2; /*!< Map of a string and a TH2D histogram, for easy 2D histogram handling. */
    static std::unordered_map<std::string, THnSparseD * > histon; /*!< Map of a string and a THnSparseD histogram, for easy nSparse handling. */
    static std::unordered_map<std::string, TNtupleD * > ttupple; /*!< Map of a string and a TNtupleD histogram, for easy Ntuple handling. */
    static std::unordered_map<std::string, TTree * > trees; /*!< Map of a string and a TTree histogram, for easy tree handling. */
    static std::unordered_map<std::string, TEfficiency * > effs; /*!< Map of a string and a TEfficiency container. */
    static std::unordered_map<std::string, TProfile * > prof; /*!< Map of a string and a TEfficiency container. */

    /*! \brief Function to create a number of 1D histograms in the histo map
     *
     * \param[in] n_histos Number of histograms that should be created with different numbers
     * \param[in] name Name of the histograms that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(int n_histos, const char* name, int nbinsx, double xlow, double xup, TString xtitle = "") {
        std::string _name=name;
        for (int j = 0; j < n_histos; j++) {
            TH1D * tmphist = new TH1D(Form("h1_%d_%s", j, _name.c_str()), xtitle, nbinsx, xlow, xup);
            tmphist->SetXTitle(xtitle);
            tmphist->Sumw2();
            histo[Form("h1_%d_%s", j, _name.c_str())] = tmphist;
        }
    }

    /*! \brief Function to create a number of 1D histograms in the histo map for a specific particle
     *
     * \param[in] n_histos Number of histograms that should be created with different numbers
     * \param[in] name Name of the histograms that should be created
     * \param[in] particle Name of the particle for which the histograms are created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(int n_histos, const char* name, const char* particle, int nbinsx, double xlow, double xup, TString xtitle = "") {
        std::string _name=name;
        for (int i = 0; i < n_histos; i++) {
            TH1D * tmphist = new TH1D(Form("h1_%d_%s_%s", i, particle, _name.c_str()), xtitle, nbinsx, xlow, xup);
            tmphist->SetXTitle(xtitle);
            tmphist->Sumw2();
            histo[Form("h1_%d_%s_%s", i, particle, _name.c_str())] = tmphist;
        }
    }

    /*! \brief Function to create a 1D histogram in the histo map without the standard histogram naming convention (h1_...)
     *
     * \param[in] name Name of the histograms that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHistoUnchangedName(const char* name, int nbinsx, double xlow, double xup, TString xtitle = "") {
        TH1D * tmphist = new TH1D(Form("%s", name), xtitle, nbinsx, xlow, xup);
        tmphist->SetXTitle(xtitle);
        tmphist->Sumw2();
        histo[Form("%s", name)] = tmphist;
    }

    /*! \brief Function to create one 1D histograms in the histo map
     *
     * \param[in] name Name of the histogram that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] ytitle Optinal title of the y-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(const char* name, int nbinsx, double xlow, double xup, TString xtitle = "", TString ytitle = "") {
        TH1D * tmphist = new TH1D(Form("h1_%s", name), xtitle, nbinsx, xlow, xup);
        tmphist->SetXTitle(xtitle);
        tmphist->SetYTitle(ytitle);
        tmphist->Sumw2();
        histo[Form("h1_%s", name)] = tmphist;
    }

    /*! \brief Function to create one 1D histogram in the histo map for a specific particle
     *
     * \param[in] name Name of the histograms that should be created
     * \param[in] particle Name of the particle for which the histograms are created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(const char* name, const char* particle, int nbinsx, double xlow, double xup, TString xtitle = "") {
        TH1D * tmphist = new TH1D(Form("h1_%s_%s", particle, name), xtitle, nbinsx, xlow, xup);
        tmphist->SetXTitle(xtitle);
        tmphist->Sumw2();
        histo[Form("h1_%s_%s", particle, name)] = tmphist;
    }

    /*! \brief Function to create one 1D histograms in the histo map with a boost name as input
     *
     * \param[in] name Name of the histogram that should be created (boost::basic_format)
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(boost::basic_format<char> name, int nbinsx, double xlow, double xup, TString xtitle = "") {
        CreateHisto(str(name).c_str(), nbinsx, xlow, xup, xtitle);
    }

    /*! \brief Function to create one 2D histograms in the histo map
     *
     * \param[in] name Name of the histogram that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] nbinsy Number of bins on the y-axis
     * \param[in] ylow Lower edge of the y-axis
     * \param[in] yup Upper edge of the y-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] ytitle Optinal title of the y-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(const char* name, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, TString xtitle = "", TString ytitle = "") {
        std::string dummy = Form("h2_%s", name);
        histo2[dummy] = new TH2D(Form("h2_%s", name), Form("h2_%s", name), nbinsx, xlow, xup, nbinsy, ylow, yup);
        histo2[dummy] -> Sumw2();
        histo2[dummy] -> GetXaxis() -> SetTitle(xtitle);
        histo2[dummy] -> GetYaxis() -> SetTitle(ytitle);
    }

    /*! \brief Function to create a number of 2D histograms in the histo map
     *
     * \param[in] n_histos Number of histograms that should be created with different numbers
     * \param[in] name Name of the histogram that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] nbinsy Number of bins on the y-axis
     * \param[in] ylow Lower edge of the y-axis
     * \param[in] yup Upper edge of the y-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] ytitle Optinal title of the y-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(int n_histos, const char* name, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, TString xtitle = "", TString ytitle = "") {
        std::string _name=name;
        for (int i = 0; i < n_histos; i++) {
            TH2D * tmphist = new TH2D(Form("h2_%d_%s", i, _name.c_str()), Form("h2_%d_%s", i, _name.c_str()), nbinsx, xlow, xup, nbinsy, ylow, yup);
            tmphist->SetXTitle(xtitle);
            tmphist->SetYTitle(ytitle);
            tmphist->Sumw2();
            histo2[Form("h2_%d_%s", i, _name.c_str())] = tmphist;
        }
    }

    /*! \brief Function to create one 2D histograms in the histo map explicitly specifying x bin borders
     *
     * \param[in] name Name of the histogram that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xbins array with lower edge of x-axis bins and upper edge of last bin
     * \param[in] nbinsy Number of bins on the y-axis
     * \param[in] ylow Lower edge of the y-axis
     * \param[in] yup Upper edge of the y-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] ytitle Optinal title of the y-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(const char* name, int nbinsx, const double * xbins,  int nbinsy, double ylow, double yup, TString xtitle = "", TString ytitle = "") {
        std::string dummy = Form("h2_%s", name);
        histo2[dummy] = new TH2D(Form("h2_%s", name), Form("h2_%s", name), nbinsx, xbins, nbinsy, ylow, yup);
        histo2[dummy] -> Sumw2();
        histo2[dummy] -> GetXaxis() -> SetTitle(xtitle);
        histo2[dummy] -> GetYaxis() -> SetTitle(ytitle);
    }

    /*! \brief Function to create one 2D histograms in the histo map explicitly specifying y bin borders
     *
     * \param[in] name Name of the histogram that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] nbinsy Number of bins on the y-axis
     * \param[in] ybins array with lower edge of y-axis bins and upper edge of last bin
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] ytitle Optinal title of the y-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(const char* name, int nbinsx, double xlow, double xup, int nbinsy, const double * ybins, TString xtitle = "", TString ytitle = "") {
        std::string dummy = Form("h2_%s", name);
        histo2[dummy] = new TH2D(Form("h2_%s", name), Form("h2_%s", name), nbinsx, xlow, xup, nbinsy, ybins);
        histo2[dummy] -> Sumw2();
        histo2[dummy] -> GetXaxis() -> SetTitle(xtitle);
        histo2[dummy] -> GetYaxis() -> SetTitle(ytitle);
    }

    /*! \brief Function to create one 2D histograms in the histo map explicitly specifying x and y bin borders
     *
     * \param[in] name Name of the histogram that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xbins array with lower edge of x-axis bins and upper edge of last bin
     * \param[in] nbinsy Number of bins on the y-axis
     * \param[in] ybins array with lower edge of y-axis bins and upper edge of last bin
     * \param[in] yup Upper edge of the y-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] ytitle Optinal title of the y-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateHisto(const char* name, int nbinsx, const double * xbins, int nbinsy, const double * ybins, TString xtitle = "", TString ytitle = "") {
        std::string dummy = Form("h2_%s", name);
        histo2[dummy] = new TH2D(Form("h2_%s", name), Form("h2_%s", name), nbinsx, xbins, nbinsy, ybins);
        histo2[dummy] -> Sumw2();
        histo2[dummy] -> GetXaxis() -> SetTitle(xtitle);
        histo2[dummy] -> GetYaxis() -> SetTitle(ytitle);
    }


    /*! \brief Function to create one NSparse in the nSparse map
     *
     * \param[in] name Name of the NSparse that should be created
     * \param[in] dimension Number of dimensions that the NSparse should have
     * \param[in] bins Array with the number of bins for each dimension
     * \param[in] xmin Array of the lower edge of the axis for each dimension
     * \param[in] xmax Array of the upper edge of the axis for each dimension
     * \param[in] axisTitle[] Array of the axis title for each dimension
     */
    SUPPRESS_NOT_USED_WARN static void CreateNSparse(const char* name, int dimension, int* bins, double* xmin, double* xmax, std::string axisTitle[]) {
        std::string dummy = Form("hn_%s", name);
        histon[dummy] = new THnSparseD(Form("hn_%s", name), Form("hn_%s", name), dimension, bins, xmin, xmax );
        histon[dummy] -> Sumw2();
        for (int i = 0 ; i < dimension; i++) {
            histon[dummy]->GetAxis(i)->SetTitle(axisTitle[i].c_str());
        }
    }

    /*! \brief Function to create one NtupleD in the NtupleD map
     *
     * \param[in] name Name of the NSparse that should be created
     * \param[in] varlist Colon sepereated list with the name of the branches that should be created
     * \param[in] bufsize Buffer size that the NtupleD should have
     */
    SUPPRESS_NOT_USED_WARN static void CreateTree(const char* name, const char* varlist, int bufsize) {
        std::string dummy = Form("tree_%s", name);
        ttupple[dummy] = new TNtupleD(Form("tree_%s", name), name, varlist, bufsize);
    }

    /*! \brief Function to create one Tree in the TTree map
     *
     * \param[in] m Map of the name and variable that should be matched to each branch
     * \param[in] name Name of the TTree that should be created
     */
    SUPPRESS_NOT_USED_WARN static void CreateTree(std::unordered_map< std::string , float > *m, const char * name) {
        trees[name] = new TTree(name, name);
        for (std::unordered_map< std::string , float >::iterator it = m->begin(); it != m->end(); it++) {
                trees[name]->Branch(it->first.c_str(), &(it->second), Form("%s/F", it->first.c_str()));
        }
    }

    /*! \brief Function to add a branch to one Tree in the TTree map
     *
     * \param[in] branch Name of the branch to be added to the tree
     * \param[in] Variable that should be added to the tree
     * \param[in] name Name of the TTree
     */
    SUPPRESS_NOT_USED_WARN static void AddBranch(std::string branch, float content, const char * name) {
        trees[name]->Branch(branch.c_str(), &content, Form("%s/F", branch.c_str()));
    }

    /*! \brief Function to add a branch to one Tree in the TTree map
     *
     * \param[in] branch Name of the branch to be added to the tree
     * \param[in] Variable that should be added to the tree
     * \param[in] name Name of the TTree
     */
    SUPPRESS_NOT_USED_WARN static void AddBranch(std::string branch, int content, const char * name) {
        trees[name]->Branch(branch.c_str(), &content, Form("%s/I", branch.c_str()));
    }

    /*! \brief Function to add a branch to one Tree in the TTree map
     *
     * \param[in] branch Name of the branch to be added to the tree
     * \param[in] Variable that should be added to the tree
     * \param[in] name Name of the TTree
     */
    SUPPRESS_NOT_USED_WARN static void AddBranch(std::string branch, Long64_t content, const char * name) {
        trees[name]->Branch(branch.c_str(), &content, Form("%s/L", branch.c_str()));
    }

    /*! \brief Function to add a branch to one Tree in the TTree map
     *
     * \param[in] branch Name of the branch to be added to the tree
     * \param[in] Variable that should be added to the tree
     * \param[in] name Name of the TTree
     */
    SUPPRESS_NOT_USED_WARN static void AddBranch(std::string branch, UInt_t content, const char * name) {
        trees[name]->Branch(branch.c_str(), &content, Form("%s/i", branch.c_str()));
    }

    /*! \brief Function to add a branch to one Tree in the TTree map
     *
     * \param[in] branch Name of the branch to be added to the tree
     * \param[in] Variable that should be added to the tree
     * \param[in] name Name of the TTree
     */
    SUPPRESS_NOT_USED_WARN static void AddBranch(std::string branch, ULong64_t content, const char * name) {
        trees[name]->Branch(branch.c_str(), &content, Form("%s/l", branch.c_str()));
    }

    /*! \brief Function to create one 1D Efficiency container in the eff map with variable binning
     *
     * \param[in] name Name of the Efficiency container that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xbins defines the binning bounds
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateEff(const char* name, int nbinsx, double* xbins, const char* xtitle = "") {
        TEfficiency * tmpeff = new TEfficiency(Form("eff_%s", name), Form("%s;%s;%s", name, xtitle, "#epsilon"), nbinsx, xbins);
        effs[Form("eff_%s", name)] = tmpeff;
    }

    /*! \brief Function to create one 1D Efficiency container in the eff map
     *
     * \param[in] name Name of the Efficiency container that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] weighted Option for weighted events (DEFAULT = false)
     */
    SUPPRESS_NOT_USED_WARN static void CreateEff(const char* name, int nbinsx, double xlow, double xup, const char* xtitle = "", bool weighted = false) {
        TEfficiency * tmpeff = new TEfficiency(Form("eff_%s", name), Form("%s;%s;%s", name, xtitle, "#epsilon"), nbinsx, xlow, xup);
        if (weighted) {
            tmpeff -> SetUseWeightedEvents();
        }
        effs[Form("eff_%s", name)] = tmpeff;
    }

    /*! \brief Function to create one 1D Efficiency container in the eff map from two histograms
     *
     * \param[in] name Name of the Efficiency container that should be created
     * \param[in] hist_all Histogram of all events
     * \param[in] hist_pass Histogram of all passed events
     * \param[in] weighted Bool if the histograms contain weighted events (DEFAULT = false)
     */
    SUPPRESS_NOT_USED_WARN static void CreateEff(const char* name, TH1D* hist_all, TH1D* hist_pass, bool weighted = false) {
        if (weighted and TEfficiency::CheckConsistency(*hist_pass,*hist_pass,"w")) {
            TEfficiency * tmpeff = new TEfficiency();
            tmpeff -> SetUseWeightedEvents();
            tmpeff -> SetTotalHistogram(*hist_all,"f");
            tmpeff -> SetPassedHistogram(*hist_pass,"f");
            tmpeff -> SetName(Form("eff_%s", name));
            effs[Form("eff_%s", name)] = tmpeff;
        } else if (not weighted) {
            TEfficiency * tmpeff = new TEfficiency(*hist_pass, *hist_all);
            tmpeff -> SetName(Form("eff_%s", name));
            effs[Form("eff_%s", name)] = tmpeff;
        }
    }

    /*! \brief Function to create one 2D Efficiency container in the eff map
     *
     * \param[in] name Name of the Efficiency container that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] nbinsy Number of bins on the y-axis
     * \param[in] ylow Lower edge of the y-axis
     * \param[in] yup Upper edge of the y-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] ytitle Optinal title of the y-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateEff(const char* name, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, const char* xtitle = "", const char* ytitle = "") {
        TEfficiency * tmpeff = new TEfficiency(Form("eff_%s", name), Form("%s;%s;%s;%s", name, xtitle, ytitle, "#epsilon"), nbinsx, xlow, xup, nbinsy, ylow, yup);
        effs[Form("eff_%s", name)] = tmpeff;
    }

    /*! \brief Function to create one 1D Profile container in the eff map
     *
     * \param[in] name Name of the Profile container that should be created
     * \param[in] nbinsx Number of bins on the x-axis
     * \param[in] xlow Lower edge of the x-axis
     * \param[in] xup Upper edge of the x-axis
     * \param[in] xtitle Optinal title of the x-axis (DEFAULT = "")
     * \param[in] ytitle Optinal title of the y-axis (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void CreateProf(const char* name, int nbinsx, double xlow, double xup, const char* xtitle = "", const char* ytitle = "") {
        TProfile * tmpprof = new TProfile(Form("prof_%s", name), Form("%s;%s;%s", name, xtitle, ytitle), nbinsx, xlow, xup);
        prof[Form("prof_%s", name)] = tmpprof;
    }

    /*! \brief Function to delete a histogram from the histo map
     *
     * \param[in] n_histos Number of histogram to be deleted
     * \param[in] name Name of the histogram to be deleted
     */
    SUPPRESS_NOT_USED_WARN static void DeleteHisto(int n_histos, const char* name){
        for(int i = 0; i < n_histos; i++){
            std::string dummy = Form("h1_%d_%s", i, name);
            delete histo[dummy];
        }
    }

    /*! \brief Function to rebin histograms in the histo map
     *
     * \param[in] n_histos Number of histograms with the same name
     * \param[in] name Name of the histograms that should be rebinned
     * \param[in] n_rebin Number of bins that the rebinned histogram should have
     * \param[in] bins Array of bin edges that the rebinned histogram should have
     */
    SUPPRESS_NOT_USED_WARN static void RebinHisto(int n_histos, const char* name, int n_rebin, double* bins) {
        for(int i = 0; i < n_histos; i++){
            std::string dummy = Form("h1_%d_%s", i, name);
            std::unordered_map<std::string, TH1D * >::iterator it = histo.find(dummy);
            if (it != histo.end()) {
                //it->second->Fill(value, weight);
            } else {
                std::cerr << "(Rebin) No hist: " << dummy << " in map " << n_histos << std::endl;
                for(auto h :histo){
                    std::cerr <<h.first<< std::endl;
                }
            }
            //std::cout<<histo.find(dummy)!=histo.end()<<std::endl;
            //histo[dummy]->SetName((dummy+"tmp").c_str());
            //std::cout<<histo[dummy]->GetNbinsX()<< " "<<n_rebin<<std::endl;
            //char* cdummy = Form("h1_%d_%s", i, name);
            histo[dummy] = (TH1D*)it->second->Rebin(n_rebin,dummy.c_str(), bins);
        }
    }

    /*! \brief Function to fill an event in a 1D histogram of the map
     *
     * This function fills one value with one weight for one event in one
     * specific histogram. The function also checks if the histogram exists
     * in the map, otherwise it will print an error message.
     * \param[in] n_histo Number of the histogram that should be filled
     * \param[in] name Name of the histogram which should be filled
     * \param[in] value Value that should be filled
     * \param[in] weight Weight of the event that should be filled
     */
    static void Fill(int n_histo, string name, double value, double weight) {
        std::unordered_map<std::string, TH1D * >::iterator it = histo.find("h1_"+to_string(n_histo)+"_"+name);
        
        //Form("h1_%d_%s", n_histo, name));
        if (it != histo.end()) {
            it->second->Fill(value, weight);
        } else {
            std::cerr << "(Fill) No hist: " << "h1_"+to_string(n_histo)+"_"+name << " in map " << n_histo <<" size is : "<<  histo.size()<<"   "<<histo.max_size()
 << std::endl;
        }
    }


    /*! \brief Function to fill an event in a 1D histogram of the map without histo number
     *
     * This function fills one value with one weight for one event in one
     * specific histogram. The function also checks if the histogram exists
     * in the map, otherwise it will print an error message.
     * \param[in] name Name of the histogram which should be filled
     * \param[in] value Value that should be filled
     * \param[in] weight Weight of the event that should be filled
     */
    SUPPRESS_NOT_USED_WARN static void Fill(const char * name, double value, double weight) {
        std::unordered_map<std::string, TH1D * >::iterator it;
        if (strcmp(name, "h_counters") == 0) {
                it = histo.find(Form("%s", name));
        } else {
                it = histo.find(Form("h1_%s", name));
        }

        if (it != histo.end()) {
            it->second->Fill(value, weight);
        } else {
            std::cerr << "(Fill) No hist: " << Form("h1_%s", name) << " in map " << std::endl;
        }
    }

    /*! \brief Function to fill an event in a 1D histogram of the map without histo number
     *
     * This function fills one value with one weight for one event in one
     * specific histogram. The function also checks if the histogram exists
     * in the map, otherwise it will print an error message.
     * \param[in] name Name of the histogram which should be filled
     * \param[in] value Value that should be filled
     * \param[in] weight Weight of the event that should be filled
     */
    SUPPRESS_NOT_USED_WARN static void FillStr(const char * name, const char * value, double weight) {
        std::unordered_map<std::string, TH1D * >::iterator it;
        if (strcmp(name, "h_counters") == 0) {
                it = histo.find(Form("%s", name));
        } else {
                it = histo.find(Form("h1_%s", name));
        }

        if (it != histo.end()) {
            it->second->Fill(value, weight);
        } else {
            std::cerr << "(Fill) No hist: " << Form("h1_%s", name) << " in map " << std::endl;
        }
    }

    /*! \brief Function to fill an event in a 1D histogram of the map without histo number
     *
     * This function fills one value with one weight for one event in one
     * specific histogram. The function also checks if the histogram exists
     * in the map, otherwise it will print an error message.
     * \param[in] name Name of the histogram which should be filled
     * \param[in] value Value that should be filled
     * \param[in] weight Weight of the event that should be filled
     */
    SUPPRESS_NOT_USED_WARN static void FillStr(int n_histo, const char * name, const char * value, double weight) {
        FillStr(Form("%d_%s",n_histo, name), value,  weight);
    }

    /*! \brief Function to fill an event in a 2D histogram of the map
     *
     * \param[in] name Name of the histogram which should be filled
     * \param[in] valuex x-value that should be filled
     * \param[in] valuey y-value that should be filled
     * \param[in] weight Weight of the event that should be filled
     */
    SUPPRESS_NOT_USED_WARN static void Fill(const char * name, double valuex, double valuey, double weight) {
        auto it =histo2.find(Form("h2_%s", name));
        if(it != histo2.end()){
            it->second->Fill(valuex, valuey, weight);
        }else{
            std::cerr << "(Fill) No h2: " << name << " in map " << std::endl;
        }
    }

    /*! \brief Function to fill an event in a 2D histogram of the map
     *
     * This function fills one value with one weight for one event in one
     * specific histogram. The function also checks if the histogram exists
     * in the map, otherwise it will print an error message.
     * \param[in] n_histo Number of the histogram that should be filled
     * \param[in] name Name of the histogram which should be filled
     * \param[in] valuex x-value that should be filled
     * \param[in] valuey y-value that should be filled
     * \param[in] weight Weight of the event that should be filled
     */
    SUPPRESS_NOT_USED_WARN static void Fill(int n_histo, const char * name, double valuex, double valuey, double weight) {
        std::unordered_map<std::string, TH2D * >::iterator it = histo2.find(Form("h2_%d_%s", n_histo, name));
        if (it != histo2.end()) {
            it->second->Fill(valuex, valuey, weight);
        } else {
            std::cerr << "(Fill) No hist: " << Form("h2_%d_%s", n_histo, name) << " in map " << n_histo << std::endl;
        }
    }

    /*! \brief Function to fill an event in a nSparse of the nSparse map
     *
     * This function fills one value with one event in one specific
     * nSparse. The function also checks if the nSparse exists
     * in the map, otherwise it will print an error message.
     * \param[in] name Name of the n which should be filled
     * \param[in] n
     * \param[in] ...
     * \todo complete the function
     */
    SUPPRESS_NOT_USED_WARN static void FillSparse(const char * name, int n, ...) {
        std::unordered_map<std::string, THnSparseD * >::iterator it = histon.find(Form("hn_%s", name));
        if (it != histon.end()) {
            std::vector <double> v;
            va_list vl;
            va_start(vl, n);
            for (int i = 0; i < n; i++) {
                v.push_back(va_arg(vl, double));
            }
            va_end(vl);
            it->second->Fill(&v[0]);
        } else {
            std::cerr << "(Fill) No hist: " << name << " in map " << std::endl;
        }
    }

    /*! \brief Function to fill an event in a NTupleD of the map
     *
     * \param[in] name Name of the NTupleD which should be filled
     * \param[in] values Array of values that should be filled
     */
    SUPPRESS_NOT_USED_WARN static void FillTree(const char * name, double* values) {
        std::string dummy = Form("tree_%s", name);
        ttupple[dummy]->Fill(values);
    }

    /*! \brief Function to fill an event in a TTree of the map
     *
     * \param[in] name Name of the TTree which should be filled
     */
    SUPPRESS_NOT_USED_WARN static void FillTree(const char * name) {
        trees[name]->Fill();
    }

    /*! \brief Function to fill an event in a efficiency container of the map
     *
     * \param[in] name Name of the histogram which should be filled
     * \param[in] valuex x-value that should be filled
     * \param[in] passed Boolean if the event passed or not
     */
    SUPPRESS_NOT_USED_WARN static void FillEff(const char * name, double valuex, bool passed) {
        auto it =effs.find(Form("eff_%s", name));
        if(it != effs.end()){
            it->second->Fill(passed, valuex);
        }else{
            std::cerr << "(Fill) No eff: " << name << " in map " << std::endl;
        }
    }

    /*! \brief Function to fill an weighted event in a efficiency container of the map
     *
     * \param[in] name Name of the histogram which should be filled
     * \param[in] valuex x-value that should be filled
     * \param[in] passed Boolean if the event passed or not
     * \param[in] weight Weight for the event to be filled
     */
    SUPPRESS_NOT_USED_WARN static void FillEff(const char * name, double valuex, bool passed, double weight) {
        auto it =effs.find(Form("eff_%s", name));
        if(it != effs.end()){
            it->second->FillWeighted(passed, weight, valuex);
        }else{
            std::cerr << "(Fill) No eff: " << name << " in map " << std::endl;
        }
    }

    /*! \brief Function to fill an event in a 2D efficiency container of the map
     *
     * \param[in] name Name of the histogram which should be filled
     * \param[in] valuex x-value that should be filled
     * \param[in] valuey y-value that should be filled
     * \param[in] passed Boolean if the event passed or not
     */
    SUPPRESS_NOT_USED_WARN static void FillEff(const char * name, double valuex, double valuey, bool passed) {
        auto it =effs.find(Form("eff_%s", name));
        if(it != effs.end()){
            it->second->Fill(passed, valuex, valuey);
        }else{
            std::cerr << "(Fill) No eff: " << name << " in map " << std::endl;
        }
    }

    /*! \brief Function to fill an event in a 1D profile container of the map
     *
     * \param[in] name Name of the histogram which should be filled
     * \param[in] valuex x-value that should be filled
     * \param[in] valuey y-value that should be filled
     * \param[in] passed Boolean if the event passed or not
     */
    SUPPRESS_NOT_USED_WARN static void Profile(const char * name, double valuex, double valuey) {
        std::string dummy = Form("prof_%s", name);
        std::unordered_map<std::string, TProfile * >::iterator it = prof.find(dummy);
        if (it != prof.end()) {
             it->second->Fill(valuex, valuey);
        } else {
            std::cerr << "(Fill) No prof: " << name << " in map " << std::endl;
        }
    }

    /*! \brief Function to write one 1D histogram of the map
     *
     * \param[in] n_histo Number of the histogram that should be written
     * \param[in] name Name of the histogram that should be written
     */
    SUPPRESS_NOT_USED_WARN static void Write(int n_histo, const char * name) {
        std::string dummy = Form("h1_%d_%s", n_histo, name);
        histo[dummy]->Write();
    }

    /*! \brief Function to write one 1D histogram of the map without number
     *
     * \param[in] name Name of the histogram that should be written
     */
    SUPPRESS_NOT_USED_WARN static void Write(const char * name) {
        std::string dummy = Form("h1_%s", name);
        histo[dummy]->Write();
    }

    /*! \brief Function to write many 1D histograms of the map
     *
     * This function writes all histograms of the map with the
     * default options, otherwise it writes all histograms that
     * contain the given string in there name.
     * \param[in] name Optional string that all histogram names that should be written contain (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void WriteAll(const char * name = "") {
        std::unordered_map<std::string, TH1D * >::iterator it;
        for (std::unordered_map<std::string, TH1D * >::iterator it = histo.begin(); it != histo.end(); ++it) {
            if (strcmp(name, "") != 0 && std::string::npos != it->first.find(name)) {
                it->second -> Write();
            } else if (strcmp(name, "") == 0) {
                it->second -> Write();
            }
        }
    }
    
    /*! \brief Function to write many 1D histograms of the map
     *
     * This function writes all histograms of the map with the
     * default options, otherwise it writes all histograms that
     * contain the given string in there name.
     * \param[in] name Optional string that all histogram names that should be written contain (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void WriteAll_NonZero(const char * name = "") {
        std::unordered_map<std::string, TH1D * >::iterator it;
        for (std::unordered_map<std::string, TH1D * >::iterator it = histo.begin(); it != histo.end(); ++it) {
            if (strcmp(name, "") != 0 && std::string::npos != it->first.find(name)) {
                if(it->second->Integral()>0){
                    it->second -> Write();
                }
            } else if (strcmp(name, "") == 0) {
                if(it->second->Integral()>0){
                    it->second -> Write();
                }
            }
        }
    }

    /*! \brief Function to write many 2D histograms of the map
     *
     * This function writes all histograms of the map with the
     * default options, otherwise it writes all histograms that
     * contain the given string in there name.
     * \param[in] name Optional string that all histogram names that should be written contain (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void WriteAll2(const char * name = "") {
        std::unordered_map<std::string, TH2D * >::iterator it;
        for (std::unordered_map<std::string, TH2D * >::iterator it = histo2.begin(); it != histo2.end(); ++it) {
            if (strcmp(name, "") != 0 && std::string::npos != it->first.find(name)) {
                it->second -> Write();
            } else if (strcmp(name, "") == 0) {
                it->second -> Write();
            }
        }
    }

    /*! Example to create a nice folder structure in your output folder
     *   //void specialAna::channel_writer(TFile* file, const char* channel) {
     *       //file1->cd();
     *       //file1->mkdir(channel);
     *       //for ( int i = 0; i < channel_stages[channel]; i++) {
     *           //char n_satge = (char)(((int)'0')+i);
     *           //file1->mkdir(TString::Format("%s/Stage_%c", channel, n_satge));
     *           //file1->cd(TString::Format("%s/Stage_%c/", channel, n_satge));
     *           //HistClass::WriteAll(TString::Format("_%s_", channel),TString::Format("%s:_%c_", channel, n_satge),TString::Format("sys"));
     *           //file1->cd();
     *           //file1->mkdir(TString::Format("%s/Stage_%c/sys", channel, n_satge));
     *           //file1->cd(TString::Format("%s/Stage_%c/sys/", channel, n_satge));
     *           //HistClass::WriteAll(TString::Format("_%s_", channel),TString::Format("_%c_:sys", n_satge));
     *       //}
     *       //file1->cd();
     *   //}
    */

    /*! \brief Function split a string at a delimiter
     *
     * This function splits a given string at a given delimineter,
     * and pushes the results in a given vector.
     * \param[in] &s String that should be split
     * \param[in] delim Delimiter where the string should be split
     * \param[in] &elems Vector in which the substrings should be pushed
    */
    SUPPRESS_NOT_USED_WARN void split(const std::string &s, char delim, std::vector<std::string> *elems) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems->push_back(item);
        }
    }

    /*! \brief Function split a string at a delimiter and return the results
     *
     * This function splits a given string at a given delimineter,
     * and returns the resulting substrings as a vector.
     * \param[in] &s String that should be split
     * \param[in] delim Delimiter where the string should be split
     * \param[out] elems Vector in which the substrings were pushed
    */
    SUPPRESS_NOT_USED_WARN std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split(s, delim, &elems);
        return elems;
    }

    /*! \brief Function to write many 1D histograms which contain specific strings of the map
     *
     * This function writes all histograms of the map with that
     * contain the given string in there name. The written histo-
     * grams also have to contain a list of strings that are sepe-
     * rated by a ':'.
     * \param[in] name String that all histogram names that should be written contain
     * \param[in] contains_i String that of names (seperated by ':') that the histogram name should contain
     */
    SUPPRESS_NOT_USED_WARN static void WriteAll(const char * name, const char * contains_i) {
        const std::string contains(contains_i);
        std::vector<std::string> i_cont = split(contains, ':');
        std::unordered_map<std::string, TH1D * >::iterator it;
        for (std::unordered_map<std::string, TH1D * >::iterator it = histo.begin(); it != histo.end(); ++it) {
            if (std::string::npos != it->first.find(name)) {
                bool do_write = false;
                for (uint i = 0; i < i_cont.size(); i++) {
                    if (std::string::npos != it->first.find(i_cont[i])) {
                        do_write = true;
                    } else {
                        do_write = false;
                        break;
                    }
                }
                if ( do_write ) {
                    it->second -> Write();
                }
            }
        }
    }

    /*! \brief Function to write many 2D histograms which contain specific strings of the map
     *
     * This function writes all histograms of the map with that
     * contain the given string in there name. The written histo-
     * grams also have to contain a list of strings that are sepe-
     * rated by a ':'.
     * \param[in] name String that all histogram names that should be written contain
     * \param[in] contains_i String that of names (seperated by ':') that the histogram name should contain
     */
    SUPPRESS_NOT_USED_WARN static void WriteAll2(const char * name, const char * contains_i) {
        const std::string contains(contains_i);
        std::vector<std::string> i_cont = split(contains, ':');
        std::unordered_map<std::string, TH2D * >::iterator it;
        for (std::unordered_map<std::string, TH2D * >::iterator it = histo2.begin(); it != histo2.end(); ++it) {
            if (std::string::npos != it->first.find(name)) {
                bool do_write = false;
                for (uint i = 0; i < i_cont.size(); i++) {
                    if (std::string::npos != it->first.find(i_cont[i])) {
                        do_write = true;
                    } else {
                        do_write = false;
                        break;
                    }
                }
                if ( do_write ) {
                    it->second -> Write();
                }
            }
        }
    }

    /*! \brief Function to write many 1D histograms which (not) contain specific strings of the map
     *
     * This function writes all histograms of the map with that
     * contain the given string in there name. The written histo-
     * grams also have to contain a list of strings that are sepe-
     * rated by a ':'. In this version also a list of strings that
     * should not be contained in the histogram name can be given.
     * \param[in] name String that all histogram names that should be written contain
     * \param[in] contains_i String that of names (seperated by ':') that the histogram name should contain
     * \param[in] vetos_i String that of names (seperated by ':') that the histogram name should not contain
     */
    SUPPRESS_NOT_USED_WARN static void WriteAll(const char * name, const char * contains_i, const char * vetos_i) {
        const std::string contains(contains_i);
        const std::string vetos(vetos_i);
        std::vector<std::string> i_cont = split(contains, ':');
        std::vector<std::string> i_veto = split(vetos, ':');
        std::unordered_map<std::string, TH1D * >::iterator it;
        for (std::unordered_map<std::string, TH1D * >::iterator it = histo.begin(); it != histo.end(); ++it) {
            if (std::string::npos != it->first.find(name)) {
                bool do_write = false;
                for (uint i = 0; i < i_cont.size(); i++) {
                    if (std::string::npos != it->first.find(i_cont[i])) {
                        do_write = true;
                    } else {
                        do_write = false;
                        break;
                    }
                }
                for (uint i = 0; i < i_veto.size(); i++) {
                    if (std::string::npos != it->first.find(i_veto[i])) {
                        do_write = false;
                        break;
                    }
                }
                if ( do_write ) {
                    it->second -> Write();
                }
            }
        }
    }

    /*! \brief Function to write many 2D histograms which (not) contain specific strings of the map
     *
     * This function writes all histograms of the map with that
     * contain the given string in there name. The written histo-
     * grams also have to contain a list of strings that are sepe-
     * rated by a ':'. In this version also a list of strings that
     * should not be contained in the histogram name can be given.
     * \param[in] name String that all histogram names that should be written contain
     * \param[in] contains_i String that of names (seperated by ':') that the histogram name should contain
     * \param[in] vetos_i String that of names (seperated by ':') that the histogram name should not contain
     */
    SUPPRESS_NOT_USED_WARN static void WriteAll2(const char * name, const char * contains_i, const char * vetos_i) {
        const std::string contains(contains_i);
        const std::string vetos(vetos_i);
        std::vector<std::string> i_cont = split(contains, ':');
        std::vector<std::string> i_veto = split(vetos, ':');
        std::unordered_map<std::string, TH2D * >::iterator it;
        for (std::unordered_map<std::string, TH2D * >::iterator it = histo2.begin(); it != histo2.end(); ++it) {
            if (std::string::npos != it->first.find(name)) {
                bool do_write = false;
                for (uint i = 0; i < i_cont.size(); i++) {
                    if (std::string::npos != it->first.find(i_cont[i])) {
                        do_write = true;
                    } else {
                        do_write = false;
                        break;
                    }
                }
                for (uint i = 0; i < i_veto.size(); i++) {
                    if (std::string::npos != it->first.find(i_veto[i])) {
                        do_write = false;
                        break;
                    }
                }
                if ( do_write ) {
                    it->second -> Write();
                }
            }
        }
    }

    /*! \brief Function to write many TTrees and TNtupleDs of the maps
     *
     * This function writes all TTrees and TNtupleDs of the maps
     * with the default options, otherwise it writes all histograms
     * that contain the given string in there name.
     * \param[in] name Optional string that all TTrees or TNtupleDs names that should be written contain (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void WriteAllTrees(const char * name = "") {
        for (std::unordered_map<std::string, TNtupleD * >::iterator it = ttupple.begin(); it != ttupple.end(); ++it) {
            if (strcmp(name, "") != 0 && std::string::npos != it->first.find(name)) {
                it->second -> Write();
            } else if (strcmp(name, "") == 0) {
                it->second -> Write();
            }
        }
        for (std::unordered_map<std::string, TTree * >::iterator it = trees.begin(); it != trees.end(); ++it) {
            if (strcmp(name, "") != 0 && std::string::npos != it->first.find(name)) {
                it->second -> Write();
            } else if (strcmp(name, "") == 0) {
                it->second -> Write();
            }
        }
    }

    /*! \brief Function to write many TNsparses of the map
     *
     * This function writes all nSparses of the map with the
     * default options, otherwise it writes all histograms that
     * contain the given string in there name.
     * \param[in] name Optional string that all nSparses names that should be written contain (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void WriteN(const char * name = "") {
        std::unordered_map<std::string, THnSparseD * >::iterator it;
        for (std::unordered_map<std::string, THnSparseD * >::iterator it = histon.begin(); it != histon.end(); ++it) {
            if (strcmp(name, "") != 0 && std::string::npos != it->first.find(name)) {
                it->second -> Write();
            } else if (strcmp(name, "") == 0) {
                it->second -> Write();
            }
        }
    }

    /*! \brief Function to write one 2D histogram of the map
     *
     * \param[in] name Name of the histogram that should be written
     */
    SUPPRESS_NOT_USED_WARN static void Write2(const char * name) {
        std::string dummy = Form("h2_%s", name);
        histo2[dummy]->Write();
    }

    /*! \brief Function to write many efficiency containers of the map
     *
     * This function writes all efficiency containers of the map with the
     * default options, otherwise it writes all efficiency containers that
     * contain the given string in there name.
     * \param[in] name Optional string that all efficiency containers names that should be written contain (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void WriteAllEff(const char * name = "") {
        std::unordered_map<std::string, TEfficiency * >::iterator it;
        for (std::unordered_map<std::string, TEfficiency * >::iterator it = effs.begin(); it != effs.end(); ++it) {
            if (strcmp(name, "") != 0 && std::string::npos != it->first.find(name)) {
                it->second -> Write();
            } else if (strcmp(name, "") == 0) {
                it->second -> Write();
            }
        }
    }

    /*! \brief Function to write many efficiency containers which contain specific strings of the map
     *
     * This function writes all efficiency containers of the map with that
     * contain the given string in there name. The written efficiencies also have to contain a list of strings that are sepe-
     * rated by a ':'.
     * \param[in] name String that all efficiency names that should be written contain
     * \param[in] contains_i String that of names (seperated by ':') that the efficiency name should contain
     */
    SUPPRESS_NOT_USED_WARN static void WriteAllEff(const char * name, const char * contains_i) {
        const std::string contains(contains_i);
        std::vector<std::string> i_cont = split(contains, ':');
        std::unordered_map<std::string, TH1D * >::iterator it;
        for (std::unordered_map<std::string, TEfficiency * >::iterator it = effs.begin(); it != effs.end(); ++it) {
            if (std::string::npos != it->first.find(name)) {
                bool do_write = false;
                for (uint i = 0; i < i_cont.size(); i++) {
                    if (std::string::npos != it->first.find(i_cont[i])) {
                        do_write = true;
                    } else {
                        do_write = false;
                        break;
                    }
                }
                if ( do_write ) {
                    it->second -> Write();
                }
            }
        }
    }

    /*! \brief Function to write many profile containers of the map
     *
     * This function writes all profile containers of the map with the
     * default options, otherwise it writes all profile containers that
     * contain the given string in there name.
     * \param[in] name Optional string that all profile containers names that should be written contain (DEFAULT = "")
     */
    SUPPRESS_NOT_USED_WARN static void WriteAllProf(const char * name = "") {
        std::unordered_map<std::string, TProfile * >::iterator it;
        for (std::unordered_map<std::string, TProfile * >::iterator it = prof.begin(); it != prof.end(); ++it) {
            if (strcmp(name, "") != 0 && std::string::npos != it->first.find(name)) {
                it->second -> Write();
            } else if (strcmp(name, "") == 0) {
                it->second -> Write();
            }
        }
    }

    /*! \brief Function to set all negative bin contents to zero for a 1D histogram
     *
     * \param[in] n_histo Number of the histogram that should be modified
     * \param[in] name Name of the histogram that should be modified
     */
    SUPPRESS_NOT_USED_WARN static void SetToZero(int n_histo, const char * name) {
        std::string dummy = Form("h1_%d_%s", n_histo, name);
        int Nbins2 = histo[dummy] -> GetNbinsX();
        for (int bb = 0; bb < Nbins2+1; bb++) {
            double binValue = histo[dummy] -> GetBinContent(bb);
            if (binValue < 0) {
                 histo[dummy] -> SetBinContent(bb, 0.);
            }
        }
    }

    /*! \brief Function to get one 1D histogram from the map without number
     *
     * \param[in] name Name of the histogram that should be returned
     * \param[out] histo Returned histogram
     */
    SUPPRESS_NOT_USED_WARN static TH1D* ReturnHist(const char * name) {
        std::string dummy = "";
        if (strcmp(name, "h_counters") == 0) {
            dummy = Form("%s", name);
        } else {
            dummy = Form("h1_%s", name);
        }
        return histo[dummy];
    }

    /*! \brief Function to give one 1D histogram from the map alphanumeric bin labels without number
     *
     * \param[in] name Name of the histogram that should get bin names
     * \param[in] n_bins of bins that should be renamed
     * \param[in] d_mydisc Array with the names that th bins should get
     */
    SUPPRESS_NOT_USED_WARN static void NameBins(const char * name, const uint n_bins, TString* d_mydisc) {
        std::string dummy = Form("h1_%s", name);
        for (uint i = 0; i < n_bins; i++) {
            histo[dummy]->GetXaxis()->SetBinLabel(i+1, d_mydisc[i]);
        }
    }

    /*! \brief Function to give one 2D histogram from the map alphanumeric bin labels without number
     *
     * \param[in] name Name of the histogram that should get bin names
     * \param[in] n_bins_x Number of x-bins that should be renamed
     * \param[in] x_bin_names Array with the names that the x-bins should get
     * \param[in] n_bins_y Number of y-bins that should be renamed
     * \param[in] y_bin_names Array with the names that the y-bins should get
     */
    SUPPRESS_NOT_USED_WARN static void NameBins(const char * name, const uint n_bins_x, TString* x_bin_names, const uint n_bins_y, TString* y_bin_names) {
        std::string dummy = Form("h2_%s", name);
        for (uint i = 0; i < n_bins_x; i++) {
            histo2[dummy]->GetXaxis()->SetBinLabel(i+1, x_bin_names[i]);
        }
        for (uint i = 0; i < n_bins_y; i++) {
            histo2[dummy]->GetYaxis()->SetBinLabel(i+1, y_bin_names[i]);
        }
    }

    /*! \brief Function to give one 1D histogram from the map alphanumeric bin labels
     *
     * \param[in] n_histo Number of the histogram that should get bin names
     * \param[in] name Name of the histogram that should get bin names
     * \param[in] n_bins of bins that should be renamed
     * \param[in] d_mydisc Array with the names that th bins should get
     */
    SUPPRESS_NOT_USED_WARN static void NameBins(int n_histo, const char * name, const uint n_bins, TString* d_mydisc) {
        for (int i = 0; i < n_histo; i++) {
            std::string dummy = Form("h1_%d_%s", i, name);
            for (uint i = 0; i < n_bins; i++) {
                histo[dummy]->GetXaxis()->SetBinLabel(i+1, d_mydisc[i]);
            }
        }
    }

    /*! \brief Function to clean up the memory usage of the HistClass
     *
     */
    SUPPRESS_NOT_USED_WARN static void CleanUp() {
        for (std::unordered_map<std::string, TEfficiency * >::iterator it = effs.begin(); it != effs.end(); ++it) {
            delete it->second;
        }
        for (std::unordered_map<std::string, TH1D * >::iterator it = histo.begin(); it != histo.end(); ++it) {
            delete it->second;
        }
        for (std::unordered_map<std::string, TH2D * >::iterator it = histo2.begin(); it != histo2.end(); ++it) {
            delete it->second;
        }
        for (std::unordered_map<std::string, THnSparseD * >::iterator it = histon.begin(); it != histon.end(); ++it) {
            delete it->second;
        }
        for (std::unordered_map<std::string, TNtupleD * >::iterator it = ttupple.begin(); it != ttupple.end(); ++it) {
            delete it->second;
        }
        for (std::unordered_map<std::string, TTree * >::iterator it = trees.begin(); it != trees.end(); ++it) {
            delete it->second;
        }
    }
}  // namespace HistClass
#endif
