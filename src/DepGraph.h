#ifndef DepGraph_h
#define DepGraph_h

#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include "Cut_enum.h"
#include <unordered_set>


using namespace std;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> mygraph;

class DepGraph {

public:

  DepGraph();
  void loadCuts(vector<CUTS>);
  void loadCuts(CUTS);
  bool isPresent(CUTS);
  unordered_set<int> getCuts();
  
private:
  void dfs(int vertex);
  mygraph g;

  unordered_set<int> neededCuts;
};

#endif
