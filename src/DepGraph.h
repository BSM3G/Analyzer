#ifndef DepGraph_h
#define DepGraph_h

#include <iostream>
#include "Cut_enum.h"
#include <unordered_set>
#include <boost/graph/adjacency_list.hpp>


//using namespace std;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> mygraph;

class DepGraph {

public:

  DepGraph();
  void loadCuts(std::vector<CUTS>);
  void loadCuts(CUTS);
  bool isPresent(CUTS);
  std::unordered_set<int> getCuts();
  
private:
  void dfs(int vertex);
  mygraph g;

  std::unordered_set<int> neededCuts;
};

#endif
