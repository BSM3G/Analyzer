#ifndef DATA_BINNER_H_
#define DATA_BINNER_H_

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <cassert>
#include <TH1.h>
#include <TH2.h>
#include <TEfficiency.h>
#include <TFile.h>

//using namespace std;

/*
DataPiece: parent class for Piece1D and Piece2D.  Use this polymorphism to make
using the Databinner a little easier

Contains the public functions for the DataBinner to use:

write_histogram(std::vector<std::string>& folders, TFiles* outfile)
Writes the histograms stored in the DataPiece into each folder of the outfile.

input: folder -- a std::vector with the names of the folders so the programs can cd into the correct folder
in the outfile
outfile -- TFile that the histograms are written out to

bin(int folder, double x, (double y), double weight)
Takes in a value and increments the correct histogram by the weight given

input: folder -- folder of histogram that will be incremented
x (y) -- value of the x axis (and y axis if 2D)
weight -- weight given to the value.

*/
class DataPiece {
protected:
  const std::string name;

public:

  bool is1D;
  DataPiece(std::string _name, int _Nfold) : name(_name), is1D(true) {
    TH1::AddDirectory(false);  
  };
  virtual ~DataPiece() {};
  virtual void write_histogram(std::vector<std::string>&, TFile*, std::string subfolder) {};
  virtual void bin(int, double, double) {};
  virtual void bin(int, double, double, double) {};
  virtual void bin(int, double, bool) {};

};


class Piece1D : public DataPiece {
private:
  const double begin, end;
  const int bins;

  std::vector<TH1D> histograms;


public:
  Piece1D(std::string, int, double, double, int);
  void write_histogram(std::vector<std::string>&, TFile*, std::string subfolder);
  void bin(int, double, double);
};


class Piece2D : public DataPiece {
private:
  const double beginx, endx, beginy, endy;
  const int binx, biny;

  std::vector<TH2D> histograms;

public:
  Piece2D(std::string, int, double, double, int, double, double, int);
  void write_histogram(std::vector<std::string>&, TFile*, std::string subfolder);
  void bin(int, double, double, double);
};


class Piece1DEff : public DataPiece {
private:
  const double begin, end;
  const int bins;
  std::vector<TEfficiency> histograms;
  bool wroteOutput;


public:
  Piece1DEff(std::string, int, double, double, int);
  void write_histogram(std::vector<std::string>&, TFile*);
  void bin(int, double, bool);
};


/*
  DataBinner: Class that is a container for the DataPieces.  All of its functions are simply to interface with individual
       DataPieces.  Uses an std::unordered_map because they have O(1) look up time (as opposed to regular std::map with O(logn)), but 
       have to include the order to write out.  Speed lost in this, but since write out happens only once, speed should
       be gained here

  AddPoint(std::string shortname, int 
  Add_Hist(std::string shortname, std::string fullname, int bin, double left, double right, int Nfolder)
       Function takes in parameters for a histogram and makes a corrisponding DataPiece for this histogram.  
       Two versions for Piece1D and Piece2D

       input: shortname -- This is the name used to identify the histogram generally.  The short name is gerenated by 
                           removing the particle identifies.  This allows the code to be general, eg if looking for the 
                           Pt of several particles, instead of labeling each "MuonPt", "TauPt", etc, they are all labelled 
			   "Pt" and the particle is infered by the group name of the particle which has the full name.
                           More specifics on the shortname and how it's generated in the github wiki page
			   https://github.com/BSM3G/Analyzer/wiki
               fullname -- The full name of the histogram.  This is what will appear as the title in the outfile
                 bin    -- The number of bins in the histogram
                left    -- The lower limit of the histogram
                right   -- The upper limit of the histogram
               Nfolder  -- The number of folders in the outfile.  Used for writing same histogram (eg MET) to different
	                   folders with different cuts
 */
class DataBinner {
public:
  DataBinner();
  DataBinner(const DataBinner&);
  DataBinner(DataBinner&&);
  DataBinner& operator=(const DataBinner&);
  DataBinner& operator=(DataBinner&&);
  ~DataBinner();

  void AddPoint(std::string,int, double, double);
  void AddPoint(std::string,int, double, double, double);
  void Add_Hist(std::string, std::string, int, double, double, int);
  void Add_Hist(std::string, std::string, int, double, double, int, double, double, int);
  void Add_Hist(std::string, int, double, double, int);
  void AddEff(std::string, int, double, bool);
  void write_histogram(TFile*, std::vector<std::string>&, std::string);
  void setSingleFill() {fillSingle = true;}

private:
  std::unordered_map<std::string, DataPiece*> datamap;
  std::vector<std::string> order;
  bool fillSingle = false;
};

#endif
