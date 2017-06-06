#include "DataBinner.h"

using namespace std;

Piece1D::Piece1D(string _name, int _bins, double _begin, double _end, int _Nfold) :  
  DataPiece(_name, _Nfold, (_bins+2)), begin(_begin), end(_end), bins(_bins), width((end - begin)/bins) {}


int Piece1D::get_bin(double y) const {
  return (int)((y-begin)/width);
}

void Piece1D::bin(int folder, double y, double weight) {
  int by = (y < begin) ? -1 : get_bin(y);
  by = (y > end) ? bins+1 : by+1;
  data[fold_width*folder + by] += weight;
}

void Piece1D::write_histogram(vector<string>& folders, TFile* outfile) {
  double entries;
  TH1D histogram(name.c_str(), name.c_str(), bins, begin, end);
  
  for(int i =0; i < (int)folders.size(); i++) {
    outfile->cd(folders.at(i).c_str());
    entries = 0;
    for(int j = 0; j < (bins+2); j++) {
      histogram.SetBinContent(j, data[i*fold_width + j]);
      entries += data[i*fold_width + j];
    }
    histogram.SetEntries((int)(entries+0.5));
    histogram.Write();

  }
}

/*------------------------------------------------------------------------------------------*/

Piece2D::Piece2D(string _name, int _binx, double _beginx, double _endx, int _biny, double _beginy, double _endy, int _Nfold) :
  DataPiece(_name, _Nfold, (_binx+2)*(_biny+2)), beginx(_beginx), endx(_endx), beginy(_beginy), endy(_endy), binx(_binx), biny(_biny),  widthx((endx - beginx)/binx), widthy((endy - beginy)/biny) {
  
  is1D = false;
}


int Piece2D::get_bin(double val, bool isX) const {
  if(isX) return (int)((val-beginx)/widthx);
  else return (int)((val-beginy)/widthy);
  return -1;
}

void Piece2D::bin(int folder, double x, double y, double weight) {
  int bx = (x < beginx) ? -1 : get_bin(x, true);
  bx = (x > endx) ? binx+1 : bx+1;
  int by = (y < beginy) ? -1 : get_bin(y, false);
  by = (y > endy) ? biny+1 : by+1;
  data[fold_width*folder + bx + by*(binx+2)] += weight;
  
}

void Piece2D::write_histogram(vector<string>& folders, TFile* outfile) {
  double entries;
  TH2D histogram(name.c_str(), name.c_str(), binx, beginx, endx, biny, beginy, endy);

  for(size_t i =0; i < folders.size(); i++) {
    
    outfile->cd(folders.at(i).c_str());
    entries = 0;
    for(int j = 0; j < (binx+2)*(biny+2); j++) {
      histogram.SetBinContent(j, data[i*fold_width + j]);
      entries += data[i*fold_width + j];
    }

    histogram.SetEntries((int)(entries+0.5));
    histogram.Write();
  }
}

/*---------------------------------------------------------------------------------------*/			      
 
DataBinner::DataBinner(){}

DataBinner::DataBinner(const DataBinner& rhs) : CR(rhs.CR) {
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

DataBinner::DataBinner(DataBinner&& rhs) : CR(rhs.CR) {
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

  if(CR) {
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

  if(CR) {
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

