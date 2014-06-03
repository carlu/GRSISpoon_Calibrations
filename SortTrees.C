//To compile:
// g++ SortTrees.C CoincEff.C Calib.C PropXtalk.C CalibTools.C Options.C Utils.C -I$GRSISYS/include --std=c++0x -o SortTrees $GRSISYS/libraries/TigFormat/libFormat.so $GRSISYS/libraries/libCalManager.so $GRSISYS/libraries/libRootIOManager.so -O0 `root-config --cflags --libs` -lTreePlayer -lSpectrum -lgsl -lgslcblas -g
//To run:
// ./Sort -f InFile1 [InFile2...] [-e (energy Calibration File)] [-w (Wave calibration file)] [-s (Source)]
// --------------------------------------------------------------------------------
// ----  Sort code for processing ROOT TTrees of TTigFragments                 ----
// ----  (as produced by GRSISPoon code from TIGRESS data)                     ----
// ----                                                                        ----
// ----  Main loop of files/trees and then events is here, sorting code        ----
// ----  included from other files.                                            ----
// ----    - C Unsworth                                                        ----
// --------------------------------------------------------------------------------

// C/C++ libraries:
#include <iostream>
#include <unordered_set>
#include <vector>
using namespace std;
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <fstream>
#include <algorithm>

// ROOT libraries:
#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TCutG.h>
#include <TTreeIndex.h>
#include <TTreePlayer.h>
#include <TChain.h>
//#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TRandom3.h>           // TRandom3 is less correlated than TRandom and almost as fast.
#include <TCanvas.h>

// GRSISpoon libraries
#include "TTigFragment.h"

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "Utils.h"



// For tracking real time
TStopwatch watch;

TApplication *App;              // Pointer to root environment for plotting etc

// Counters
int FileCount = 0;
int TreeEventCount = 0;
int ChainEventCount = 0;
int TreeFragCount = 0;
int ChainFragCount = 0;
int EmptyEventCount = 0;

// Storing alternate Calibration
vector < string > EnCalibNames;
vector < vector < float >>EnCalibValues;
// Waveform calibration
vector < string > WaveCalibNames;
vector < vector < float >>WaveCalibValues;

// Functions
//int LoadDefaultSettings();
//int ReadCommandLineSettings(int argc, char **argv);
//void PrintHelp();

void SortTree(const char *fn);
void IncSpectra();
int ReadCalibrationFile(std::string filename, vector < string > *EnCalibNames,
                        vector < vector < float >>*EnCalibValues);

void CoincEff(std::vector < TTigFragment > &ev);
int InitCoincEff();
void FinalCoincEff();

void Calib(std::vector < TTigFragment > &ev);
int InitCalib();
void FinalCalib();

int CalibSpectra(std::string filename);

void PropXtalk(std::vector < TTigFragment > &ev);
int InitPropXtalk();
void FinalPropXtalk();

void GeTiming(std::vector < TTigFragment > &ev);
int InitGeTiming();
void FinalGeTiming();

int main(int argc, char **argv)
{
   // Variables, Constants, etc
   int i;
   unsigned int ChainEvent, TreeEvent;
   unsigned int FileNum;

   // Set default and read custom options
   LoadDefaultSettings();
   if (ReadCommandLineSettings(argc, argv) < 0) {
      cout << "Failed to configure the run - exiting!" << endl;
      return -1;
   }

   // Set options for histo stats
   gStyle->SetOptStat("iouRMen");

   // create root environment for interacting with plots etc
   App = new TApplication("Output", 0, NULL);

   // Timing
   TStopwatch StopWatch;
   StopWatch.Start();

   // Load any alternate calibration information 
   if (Config.HaveAltEnergyCalibration) {
      int NumCal;
      NumCal = ReadCalibrationFile(Config.EnergyCalibrationFile, &EnCalibNames, &EnCalibValues);
      if (Config.PrintBasic) {
         cout << "Alternate energy calibratrion values: " << endl;
         for (i = 0; i < NumCal; i++) {
            cout << i;
            cout << ": " << EnCalibNames.at(i);
            cout << " " << EnCalibValues.at(i)[0] << " " << EnCalibValues[i][1] << " " << EnCalibValues[i][2] << endl;
         }
      }
   }
   // Load gain coefficients for waveforms
   if (Config.HaveWaveCalibration) {
      int NumCal;
      NumCal = ReadCalibrationFile(Config.WaveCalibrationFile, &WaveCalibNames, &WaveCalibValues);
      if (Config.PrintBasic) {
         cout << "Wave energy calibratrion values: " << endl;
         for (i = 0; i < NumCal; i++) {
            cout << i;
            cout << ": " << WaveCalibNames.at(i);
            cout << " " << WaveCalibValues.at(i)[0] << " " << WaveCalibValues[i][1] << " " << WaveCalibValues[i][2] <<
                endl;
         }
      }
   }
   // Initialise spectra   
   if (Config.RunEfficiency) {
      if (Config.PrintBasic) {
         cout << "Initialising Efficiency Spectra..." << endl;
      }
      if (InitCoincEff() != 0) {
         cout << "InitCoincEff Failed!" << endl;
         return 1;
      }
   }
   if (Config.RunCalibration) {
      if (Config.PrintBasic) {
         cout << "Initialising Calibration Spectra..." << endl;
      }
      if (InitCalib() != 0) {
         cout << "InitCalib Failed!" << endl;
         return 1;
      }
   }
   if (Config.RunPropCrosstalk) {
      if (Config.PrintBasic) {
         cout << "Initialising Cross-talk Spectra..." << endl;
      }
      if (InitPropXtalk() != 0) {
         cout << "InitPropXtalk Failed!" << endl;
         return 1;
      }
   }
   if (Config.RunGeTiming) {
      if (Config.PrintBasic) {
         cout << "Initialising Ge Timing Spectra..." << endl;
      }
      if (InitGeTiming() != 0) {
         cout << "InitGeTiming Failed!" << endl;
         return 1;
      }
   }
   

   TChain *Chain = new TChain("FragmentTree");

   // Setup input file list                                        
   for (FileNum = 0; FileNum < Config.files.size(); FileNum++) {
      Chain->Add(Config.files.at(FileNum).c_str());
      FileCount++;
   }

   //TTigFragment *pFrag = 0;  
   // changed above line to one below trying to fix memory leak when looping chain.  
   // It didn't work but root website suggests doing it with"new" so I will stick with it for now.
   TTigFragment *pFrag = new TTigFragment();

   std::vector < TTigFragment > evFrags;

   int nTrees = Chain->GetNtrees();
   unsigned int NumChainEntries = Chain->GetEntries();
   unsigned int NumChainEvents = (int) Chain->GetMaximum("TriggerId");  // This doesn't work, TrigID reset for each tree on chain.

   if (Config.PrintBasic) {
      cout << "Chain Entries (frags)               : " << NumChainEntries << endl;
      cout << "Chain Events (Max \"TriggerId\"))   : " << NumChainEvents << endl;
   }

   int TreeNum = -1;
   int LastTreeNum = -1;
   unsigned int NumTreeEntries = 0;
   unsigned int NumTreeEvents = 0;
   int FirstTreeEvent = -1;

   for (ChainEvent = 0; ChainEvent < NumChainEntries; ChainEvent++) {

      Chain->LoadTree(ChainEvent);
      TreeNum = Chain->GetTreeNumber();

      if (TreeNum != LastTreeNum) {
         if (Config.PrintBasic) {
            cout << "Switching to TreeNum " << TreeNum + 1 << " from " << LastTreeNum +
                1 << " at chain entry " << ChainEvent << endl;
         }
         LastTreeNum = TreeNum;
      } else {
         continue;              // This works in conjuncion with the line below "i += (NumTreeEntries - 10);"
         // if i still falls in the same tree, then it is incremented without sorting until the next tree
      }

      TTree *Tree = Chain->GetTree();

      if (!Tree->GetTreeIndex()) {
         if (Config.PrintBasic) {
            printf("\nTree Index not found, Building index...");
         }
         fflush(stdout);
         Tree->BuildIndex("TriggerId", "FragmentId");
         if (Config.PrintBasic) {
            printf("  Done\n");
         }
         fflush(stdout);
      }

      TBranch *Branch = Tree->GetBranch("TTigFragment");

      Branch->SetAddress(&pFrag);
      Tree->SetMaxVirtualSize(Config.ROOT_MaxVirtSize);
      Branch->LoadBaskets();

      NumTreeEntries = Tree->GetEntries();
      FirstTreeEvent = (int) Tree->GetMinimum("TriggerId");
      NumTreeEvents = ((int) Tree->GetMaximum("TriggerId")) - FirstTreeEvent;

      TreeFragCount = 0;
      TreeEventCount = 0;
      for (TreeEvent = 0; TreeEvent < (NumTreeEvents + FirstTreeEvent); TreeEvent++) {
         //for (int TreeEvent = FirstTreeEvent; TreeEvent < NumTreeEvents; TreeEvent++) {   
         evFrags.clear();
         int FragNum = 1;

         while (Tree->GetEntryWithIndex(TreeEvent, FragNum++) != -1) {
            evFrags.push_back(*pFrag);
            TreeFragCount++;
            ChainFragCount++;
            if (DEBUG_TREE_LOOP) {
               cout << "FragNum: " << FragNum << " j: " << TreeEvent;
            }
         }
         if (DEBUG_TREE_LOOP) {
            cout << endl;
         }

         if (FragNum < 3) {     // Init to 1, incremented on first GetEntryWithIndex, so should be at least 3 if any frags found      
            EmptyEventCount++;
            continue;
         } else {
            TreeEventCount++;
            ChainEventCount++;
         }

         // do something with the evFrag vector (contains a built event)... ProcessEvent(evFrags); 
         if (Config.EventLimit > 0 && ChainEventCount >= Config.EventLimit) {
            if (Config.PrintBasic) {
               cout << "Maximum number of events (" << Config.EventLimit << ") reached.  Terminating..." << endl;
            }
            break;
         }
         if (Config.RunEfficiency) {
            CoincEff(evFrags);
         }                      //passing vector of built events.
         if (Config.RunCalibration) {
            Calib(evFrags);
         }
         if (Config.RunPropCrosstalk) {
            PropXtalk(evFrags);
         }
         if(Config.RunGeTiming) {
            GeTiming(evFrags);
         }
         
         if (DEBUG_TREE_LOOP) {
            cout << "ev Num =  " << TreeEvent << ", ev.size() = " << evFrags.size() << endl;
         }
         // Print info to stdout
         if (Config.PrintBasic && (TreeEvent % PRINT_FREQ) == 0) {
            cout << "----------------------------------------------------------" << endl;
            cout << "  Tree  " << TreeNum + 1 << " / " << nTrees << "  (Frag ";
            cout << ChainFragCount << " / " << NumChainEntries << "):" << endl;
            cout << "\tEvent " << TreeEventCount << " / " << NumTreeEvents << endl;
            cout << "\tFrag  " << TreeFragCount << " / " << NumTreeEntries << endl;
            cout << "  Time: " << StopWatch.RealTime() << " seconds." << endl;
            StopWatch.Continue();
            cout << "----------------------------------------------------------" << endl;
         }

      }

      if (Config.EventLimit > 0 && ChainEventCount >= Config.EventLimit) {
         break;
      }

      Branch->DropBaskets("all");       // Clear cache before next tree    
      //i += (nEntries - 10);
      TreeNum += (NumTreeEntries - 10); // This skips i to almost the end of the tree, any remaining entries
      // on this tree will be skipped by "if(TreeNum != LastTreeNum" condition above
   }

   // Now finalise sorts and write spectra to files....  
   if (Config.RunEfficiency) {
      FinalCoincEff();
   }
   if (Config.RunCalibration) {
      FinalCalib();
   }
   if (Config.RunPropCrosstalk) {
      FinalPropXtalk();
   }
   if (Config.RunGeTiming) {
      FinalGeTiming();
   }

   return 0;
}


void IncSpectra()
{                               // For testing only
   //hTestSpectrum->SetBinContent(1,1000);
}
