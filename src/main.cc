#include "Analyzer.h"
#include <csignal>

bool do_break;
void KeyboardInterrupt_endJob(int signum) {
    do_break = true;
}

void usage() {
  cout << "./Analyzer infile.root outfile.root\n";
  cout << "or\n";
  cout << "./Analyzer -out outfile.root -in infile.root infile.root infile.root ...\n";
  cout << "or\n";
  cout << "./Analyzer -out outfile.root -in infile.root\n";
  cout << "Available options are:\n";
  cout << "-CR: to run over the control regions (not the usual output)\n";
  cout << "-C: use a different config folder than the default 'PartDet'\n";
  cout << "-t: run over 100 events\n";
  cout << "\n";

  exit(EXIT_FAILURE);
}

void parseCommandLine(int argc, char *argv[], vector<string> &inputnames, string &outputname, bool &setCR, bool &testRun, string &configFolder) {
  if(argc < 3) {
    cout << endl;
    cout << "You have entered too little arguments, please type:\n";
    usage();
  }
  for (int arg=1; arg<argc; arg++) {
    //// extra arg++ are there to move past flags
    if (strcmp(argv[arg], "-CR") == 0) {
      setCR = true;
      continue;
    }else if (strcmp(argv[arg], "-t") == 0) {
      testRun = true;
      continue;
    }else if (strcmp(argv[arg], "-C") == 0) {
      configFolder=argv[arg+1];
      cout << "Analyser: ConfigFolder " << configFolder << endl;
      arg++;
      continue;
    }else if (strcmp(argv[arg], "-in") == 0) {
      arg++;
      while( arg<argc and (argv[arg][0] != '-')){
        inputnames.push_back(argv[arg]);
        cout << "Analyser: Inputfilelist " << inputnames.back() << endl;
        arg++;
      }
      arg--; /// to counteract arg++ that is in the for loop
      continue;
    }else if (strcmp(argv[arg], "-out") == 0) {
      outputname=argv[arg+1];
      cout << "Analyser: Outputfile " << outputname << endl;
      arg++;
      continue;
    } else if(argv[arg][0] == '-') {
      cout << endl;
      cout << "You entered an option that doesn't exist.  Please use one of the options:" << endl;
      usage();
    }else if(inputnames.size()==0){
      inputnames.push_back(argv[arg]);
      cout << "Analyser: Inputfilelist " << inputnames.back() << endl;
    }else if(outputname==""){
      outputname = argv[arg];
      cout << "Analyser: Outputfile " << outputname << endl;
    }
  }

  if(inputnames.size() == 0) {
    cout << endl;
    cout << "No input files given!  Please type:" << endl;
    usage();
  } else if(outputname == "") {
    cout << endl;
    cout << "No output file given!  Please type:" << endl;
    usage();
  }


  //for( auto file: inputnames) {
    //ifstream ifile(file);
    //if ( !ifile && file.find("root://") == string::npos && file.find("root\\://") == string::npos) {
      //std::cout << "The file '" << inputnames.back() << "' doesn't exist" << std::endl;
      //exit(EXIT_FAILURE);
    //}
  //}
  return;
}

int main (int argc, char* argv[]) {

  bool setCR = false;
  bool testRun = false;
  do_break =false;

  string outputname;
  string configFolder="PartDet";
  vector<string> inputnames;


  //get the command line options in a nice loop
  parseCommandLine(argc, argv, inputnames, outputname, setCR, testRun, configFolder);


  //setup the analyser
  Analyzer testing(inputnames, outputname, setCR, configFolder);

  //catch ctrl+c and just exit the loop
  //this way we still have the output
  signal(SIGINT,KeyboardInterrupt_endJob);

  size_t Nentries=testing.nentries;
  if(testRun){
    Nentries=100;
    testing.nentries=100;
  }
  //main event loop
  for(size_t i=0; i < Nentries; i++) {
    testing.clear_values();
    testing.preprocess(i);
    testing.fill_histogram();
    //this will be set if ctrl+c is pressed
    if(do_break){
      testing.nentries=i+1;
      break;
    }
  }
  testing.printCuts();
  return 0;
}
