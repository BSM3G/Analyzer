#ifndef analysis_h
#define analysis_h

struct Analyzer;
#include "src/Analyzer.h"

class SpechialAnalysis {
public:
  SpechialAnalysis(Analyzer* _a);
  void init();

  void begin_run();
  void analyze();
  void end_run();

private:
  Analyzer* a;
  const std::string particles[4] = {"Ele", "Muon", "Tau", "MET"};
  const std::string particleSymbols[4]       = {"e", "#mu", "#tau", "E_{T}^{miss}"};
};

#endif
