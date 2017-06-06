#ifndef CRTest_h
#define CRTest_h

class Analyzer;
#include "Analyzer.h"



class CRTester {
 public:
  const FillVals* info;
  const string variable;
  const double cutVal;
  const string partName;

  CRTester(FillVals* _info, string var, double val, string _name);
  bool test(Analyzer* analyzer);
  bool partPassBoth(Analyzer*);
};


#endif
