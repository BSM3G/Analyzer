#ifndef CRTest_h
#define CRTest_h

class Analyzer;
#include "Analyzer.h"



class CRTester {
 public:
  const FillVals* info;
  const std::string variable;
  const double cutVal;
  const std::string partName;

  CRTester(FillVals* _info, std::string var, double val, std::string _name);
  bool test(Analyzer* analyzer);
  bool partPassBoth(Analyzer*);
};


#endif
