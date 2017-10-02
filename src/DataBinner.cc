#include "DataBinner.h"

using namespace std;

Piece1D::Piece1D(string _name, int _bins, double _begin, double _end, int _Nfold) :
DataPiece(_name, _Nfold), begin(_begin), end(_end), bins(_bins) {
  for(int i = 0; i < _Nfold; i++) {
    //string hname = name + to_string(i);
    string hname = name;
    TH1D tmp((hname).c_str(), name.c_str(), bins, begin, end);
    tmp.Sumw2();
    histograms.push_back(tmp);
  }
}

void Piece1D::bin(int folder, double y, double weight) {
  histograms.at(folder).Fill(y,weight);
}

void Piece1D::write_histogram(vector<string>& folders, TFile* outfile) {
  for(int i =0; i < (int)folders.size(); i++) {
    outfile->cd(folders.at(i).c_str());
    histograms.at(i).Write();
  }
}

/*------------------------------------------------------------------------------------------*/

Piece2D::Piece2D(string _name, int _binx, double _beginx, double _endx, int _biny, double _beginy, double _endy, int _Nfold) :
DataPiece(_name, _Nfold), beginx(_beginx), endx(_endx), beginy(_beginy), endy(_endy), binx(_binx), biny(_biny) {

  is1D = false;


  for(int i = 0; i < _Nfold; ++i) {
    //string hname = name + to_string(i);
    string hname = name;
    TH2D tmp((hname).c_str(), name.c_str(), binx, beginx, endx, biny, beginy, endy);
    tmp.Sumw2();
    histograms.push_back(tmp);
  }
}

void Piece2D::bin(int folder, double x, double y, double weight) {
  histograms.at(folder).Fill(x,y,weight);

}

void Piece2D::write_histogram(vector<string>& folders, TFile* outfile) {
  for(size_t i =0; i < folders.size(); i++) {
    outfile->cd(folders.at(i).c_str());
    histograms.at(i).Write();
  }
}

/*---------------------------------------------------------------------------------------*/

DataBinner::DataBinner(){}

DataBinner::DataBinner(const DataBinner& rhs) : fillSingle(rhs.fillSingle) {
  cout << "copied" << endl;
  order = rhs.order;

  for(auto it: rhs.datamap) {
    if(it.second->is1D) {
      datamap[it.first] = new Piece1D(*static_cast<Piece1D*>(it.second));
    } else {
      datamap[it.first] = new Piece2D(*static_cast<Piece2D*>(it.second));
    }
  }

}

DataBinner::DataBinner(DataBinner&& rhs) : fillSingle(rhs.fillSingle) {
  cout << "moved" << endl;
  for(auto it: datamap) {
    if(it.second != nullptr) {
      delete it.second;
      it.second = nullptr;
    }
  }

  order = rhs.order;
  datamap.swap(rhs.datamap);

  rhs.datamap.clear();
}


DataBinner::~DataBinner() {
  for(auto it: datamap) {
    if( it.second != nullptr) {
      delete it.second;
      it.second = nullptr;
    }
  }
}

void DataBinner::Add_Hist(string shortname, string fullname, int bin, double left, double right, int Nfolder) {
  datamap[shortname] = new Piece1D(fullname, bin, left, right, Nfolder);
  order.push_back(shortname);
}

void DataBinner::Add_Hist(string shortname, string fullname, int binx, double leftx, double rightx, int biny, double lefty, double righty, int Nfolder) {
  datamap[shortname] = new Piece2D(fullname, binx, leftx, rightx, biny, lefty, righty, Nfolder);
  order.push_back(shortname);
}


void DataBinner::AddPoint(string name, int maxfolder, double value, double weight) {
  if(datamap.count(name) == 0)  return;

  if(fillSingle) {
    if(maxfolder < 0) return;
    datamap.at(name)->bin(maxfolder,value, weight);
  } else {

    for(int i=0; i < maxfolder; i++) {
      datamap.at(name)->bin(i,value, weight);
    }
  }
}

void DataBinner::AddPoint(string name, int maxfolder, double valuex, double valuey, double weight) {
  if(datamap.count(name) == 0) return;

  if(fillSingle) {
    if(maxfolder < 0) return;
    datamap.at(name)->bin(maxfolder,valuex, valuey, weight);
  } else {
    for(int i=0; i < maxfolder; i++) {
      datamap.at(name)->bin(i,valuex, valuey, weight);
    }
  }
}

void DataBinner::write_histogram(TFile* outfile, vector<string>& folders) {
  for(vector<string>::iterator it = order.begin(); it != order.end(); it++) {
    datamap.at(*it)->write_histogram(folders, outfile);
  }
}
