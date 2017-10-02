# Setting Up

If you are setting up the Analyzer, click on the link for respective version

- [CMSSW_7_4_x](https://github.com/dteague/Analyzer/tree/TNT74x)
- [CMSSW_8_0_x](https://github.com/dteague/Analyzer/tree/TNT80x)

# Changes

## Version 1: Initial version, same as old BSM3G code

## Version 2: _COMING SOON_ Make MET its own cut and add SVFit
MET cuts are now in the file PartDet/Run_info.in.  SVFit cuts are not necessary for it to run.


# FAQ

- Q: [The program crashes with a SegFault](https://github.com/dteague/Analyzer#a-segfault)
- Q: [The program crashes with SegFault and Error in TTree::SetBranchStatus](https://github.com/dteague/Analyzer#segfault-with-tbranch-error)
- Q: [Why isn't my MET cut working?](#met-cut)
- Q: [How do I set up folders?](https://github.com/dteague/Analyzer#folders)
- Q: [How do I control which histograms make it into my root file?](https://github.com/dteague/Analyzer#histogram-management)
- Q: [How do I add a new histogram?](https://github.com/dteague/Analyzer#new-histograms)
- Q: [What is SVFit?  How do I use it?](https://github.com/dteague/Analyzer#svfit)

### A: SegFault

If you get a setfault, it can mean one of two things.  If the error looks like:

```
$ ./Analyzer OutTree.root test.root
setup start
TOTAL EVENTS: ###
setup complete

 *** Break *** segmentation violation

===========================================================
There was a crash.
This is the entire stack trace of all threads:
===========================================================
...
...
###  0x00007fc6e0d9d3f7 in std::__throw_out_of_range (__s=__s entry=0x47dd90 "_Map_base::at") at ../../../../../libstdc++-v3/src/c++11/functexcept.cc:90
...
```
This is a map out of bound error.  This means one of your values is not named correctly or is being parsed as the wrong values.  To check which values, look at the top of the stack.  I should look something like this:
```
The lines below might hint at the cause of the crash.
If they do not help you then please submit a bug report at
http://root.cern.ch/bugs. Please post the ENTIRE stack trace
from above as an attachment in addition to anything else
that might help us fixing this issue.
===========================================================
...
...
#12 Analyzer::getGoodRecoLeptons (this=this
entry=0x7fff69207790, lep=..., ePos=ePos
entry=CUTS::eRTau1, eGenPos=eGenPos
entry=CUTS::eGTau, stats=...) at src/Analyzer.cc:561
#13 0x0000000000463110 in Analyzer::preprocess (this=0x7fff69207790, event=0) at src/Analyzer.cc:130
#14 0x000000000041d4ac in main ()
===========================================================
```
In this example, we can see in line #12, the function called getGoodRecoLeptons, and based on the CUTS value sent in (eRTau1), we can see that the RecoTau1 has a value that is wrong.  Now we have to just go into PartDet/Tau_info.in and look under Tau1 to find the error.  To Help with this, one can look through the function in src/Analyzer.cc or look at a template info file such at in [this repository](https://github.com/dteague/Analyzer/tree/master/PartDet)

### SegFault with TBranch Error

If the Error looks like:
```
$ ./Analyzer TNT.root test.root
setup start
TOTAL EVENTS: 493
Error in <TTree::SetBranchStatus>: unknown branch -> Tau_byTightIsolationMVArun2v1DBnewDMwLT
Error in <TTree::SetBranchAddress>: unknown branch -> Tau_byTightIsolationMVArun2v1DBnewDMwLT
setup complete

 *** Break *** segmentation violation

===========================================================
There was a crash.
This is the entire stack trace of all threads:
===========================================================
...
...
...
```
The error is being thrown by ROOT because some of the Branches haven't been set correctly.  This can happen because the nTuples for 74x and 80x have different names, or because the name is simply mispelled.  ROOT tells you which branch has been set wrong so just go into PartDet/Tau_info.in and change the name there.  To find the names of the branches, simply list them by typing
```
cat NOTES
```
### Met Cut

The MET has been changed how it's implimented in the old code.  To make sure your MET is working properly, make sure you MET cut information is in the file ```PartDet/Run_info.in```.  The python script ```moveMET.py``` should do this for you if you have any doubts.

You must also tell the program to cut on the MET if required.  This is implimented in much the same was the multiplicity cuts in the file ```PartDet/Cuts.in```.  Simply put
```
METCut          1   -1
```
in the Cuts.in file to make sure the file removes events that don't pass MET cuts (or MHT and HT cuts for that matter).  Because of the change, the MET cut can be put in any order in relation to the other cuts so you can see how MET effects cut flow efficiency

### Folders

Folders in the program are made when reading PartDet/Cuts.in.  By default, the program will always make the last significant cut (range is not [0,-1]) into a folder.  To add folders, simply put ```***``` before the cut without any space.

e.g.
```
NRecoMuon1               0  -1
NRecoTau1                2   2
***NDiTauCombinations    1   0
NSusyCombinations        1  -1
NDiJetCombinations       0  -1
```
In this example, there is a cut on Tau1, DiTaus, and a VBF cut.  The folders created are NDiTauCominations and NSusyCombinations (last significant cut).

The order of the cuts can also be rearranged if one wants to see cut flow in a different way.

### Histogram Management

All of the Histograms are stored in PartDet/Hist_info.in.  On each line, the details of the histogram are:
```
<NAME>  <BINS>  <MIN>  <MAX>   // OR
<NAME 2D>  <XBINS> <XMIN>  <XMAX>  <YBINS>  <YMIN>  <YMAX>
```
Since the histogram information is read at the beginning of each run, the binning and domain of the histogram can be changed to fit the analysis.

As with all of the info files, the file supports C and python style line commenting (// and #).  This means, to remove a specific histogram, simply comment it out

Many of the histograms can be grouped together, so to facilitate the removal process, blocks of similar histograms are grouped under a heading that starts with the keywork "Fill."  To remove the block, set the Fill heading to 0 or false.  Since the heading won't be seen by the program, the calculates done by the block won't be done either, so marginal speed gains will be made by program (less 100th of total time, so not signicant)

### New Histograms

Adding a new histogram is fairly easy since the program is dynamic enough to hold most changes.  Two main things need to be done.

1. The histogram and information needs to be put into PartDet/Hist_info.in.  This includes name, bins, min, max as well as which heading the histogram will be stored under.  This can follow the template of the other histograms, so this is relatively easy
2. The histogram needs to be filled with the right values.  The filling of the histograms is done in the method ```fillFolder```.  In this method, there are several if blocks for the different headings.  Go to the appropriate heading (or make one if a new one was made in the Hist_info.in file), and write the following command to write to the histogram:
