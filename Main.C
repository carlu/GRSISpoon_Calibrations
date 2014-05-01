//To compile:
// g++ Main.C CoincEff.C Calib.C PropXtalk.C CalibTools.C CalibSpectra.C -I$GRSISYS/include --std=c++0x -o Sort $GRSISYS/libraries/TigFormat/libFormat.so $GRSISYS/libraries/libCalManager.so $GRSISYS/libraries/libRootIOManager.so -O0 `root-config --cflags --libs` -lTreePlayer -lSpectrum -lgsl -lgslcblas -g
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


/*

Medium Term:
- Add folloowing codes for sorting:
      - CoincEff.C
      - Calib.C
      - PropXtalk.C
      - Waves.C
      - DiffXtalk.C
*/

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
#include <TRandom3.h>           // TRandom3 is less correlated than TRandom and almost as fast.  ls
#include <TCanvas.h>

// GRSISpoon libraries
#include "TTigFragment.h"

// My libraries
#include "Main.h"
#include "Options.h"


// For tracking real time
TStopwatch watch;

TApplication *App;              // Pointer to root environment for plotting etc

// canvases for calibration.  To be shared between Calib() and CalibOffline() functions
TCanvas *cCalib1, *cCalib1a, *cCalib2, *cWave1, *ctemp;

// Rand for gain matching
static TRandom3 rand1;

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
int LoadDefaultSettings();
int ReadCommandLineSettings(int argc, char **argv);
void PrintHelp();

void SortTree(const char *fn);
void IncSpectra();
int ReadCalibrationFile(std::string filename, vector < string > *EnCalibNames,
                        vector < vector < float >>*EnCalibValues);

void CoincEff(std::vector < TTigFragment > &ev);
void InitCoincEff();
void FinalCoincEff();

void Calib(std::vector < TTigFragment > &ev);
void InitCalib();
void FinalCalib();

int CalibSpectra(std::string filename);

void PropXtalk(std::vector < TTigFragment > &ev);
void InitPropXtalk();
void FinalPropXtalk();

int main(int argc, char **argv)
{
   // Variables, Constants, etc
   int i, j;
   unsigned int ChainEvent, TreeEvent;


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

   // Offline calibration here
   // need to initialise TCanvas's for Calib and CalibOffline here:
   if (Config.RunSpecCal == 1 || Config.RunCalibration == 1) {
      // Condition here to test if plotting will be used

      // Initialise TCanvas's
      cCalib1 = new TCanvas("cCalib1", "Fit", 800, 600);        // Canvas for spectrum plots

      cCalib1a = new TCanvas("cCalib1a", "Calibration", 800, 600);      // Canvas for spectrum plots
      cCalib1a->Divide(1, 2);

      cCalib2 = new TCanvas("cCalib2", "Calibration Summary", 800, 600);        // Canvas for gain plots and histograms
      cCalib2->Divide(2, 3);
      cCalib2->Update();

      cCalib1->cd();
   }
   //    - at this stage, all setup required for offline calibration should be complete.
   //    - check if this is what we want to do.
   //    - check the file list for a suitable set of spectra
   //    - call the offline calibration function
   //    - return from here rather than running the TChain stuff below.
   if (Config.RunSpecCal == 1) {
      for (i = 0; i < Config.files.size(); i++) {
         if (Config.PrintBasic) {
            cout << "Attempting offline calibration on histograms in file: " << Config.files.at(i) << endl;
         }
         if (CalibSpectra(Config.files.at(i)) >= 0) {   // return after one succesful file so outputs are not overwritten.
            return 0;
         }
      }
   }
   // ---- End of histogram calibration -----

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
         //cout << WaveCalibNames.size() << endl << WaveCalibValues.size() << endl;
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
      InitCoincEff();
   }
   if (Config.RunCalibration) {
      if (Config.PrintBasic) {
         cout << "Initialising Calibration Spectra..." << endl;
      }
      InitCalib();
   }
   if (Config.RunPropCrosstalk) {
      if (Config.PrintBasic) {
         cout << "Initialising Cross-talk Spectra..." << endl;
      }
      InitPropXtalk();
   }

   TChain *Chain = new TChain("FragmentTree");

   // Setup input file list                                        
   for (i = 0; i < Config.files.size(); i++) {
      Chain->Add(Config.files.at(i).c_str());
      FileCount++;
   }

   //TTigFragment *pFrag = 0;  
   // changed above line to one below trying to fix memory leak when looping chain.  
   // It didn't work but root website suggests doing it with"new" so I will stick with it for now.
   TTigFragment *pFrag = new TTigFragment();  
   
   std::vector < TTigFragment > evFrags;

   int nTrees = Chain->GetNtrees();
   int NumChainEntries = Chain->GetEntries();
   int NumChainEvents = (int) Chain->GetMaximum("TriggerId");   // This doesn't work, TrigID reset for each tree on chain.

   if (Config.PrintBasic) {
      cout << "Chain Entries (frags) : " << NumChainEntries << endl;
      cout << "Chain Events          : " << NumChainEvents << endl;
   }

   int TreeNum = -1;
   int LastTreeNum = -1;
   int NumTreeEntries = -1;
   int NumTreeEvents = -1;
   int FirstTreeEvent = -1;
   int ProcessedFragments = 0;

   for (ChainEvent = 0; ChainEvent < NumChainEntries; ChainEvent++) {

      //cout << "HERE!" << endl;

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

      //cout << "HERE!!" << endl;

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

      //cout << "HERE!!!" << endl;

      Branch->SetAddress(&pFrag);
      Tree->SetMaxVirtualSize(Config.ROOT_MaxVirtSize);
      Branch->LoadBaskets();

      NumTreeEntries = Tree->GetEntries();
      FirstTreeEvent = (int) Tree->GetMinimum("TriggerId");
      NumTreeEvents = ((int) Tree->GetMaximum("TriggerId") ) - FirstTreeEvent;
      

      // cout << "HERE!!!! " << NumChainEntries << " " << NumTreeEntries << " " << NumTreeEvents << endl;
      TreeFragCount = 0;
      TreeEventCount = 0;
      for (int TreeEvent = 0; TreeEvent < (NumTreeEvents+FirstTreeEvent); TreeEvent++) {
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
         //cout << "HERE!!!!! " << endl;

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
      i += (NumTreeEntries - 10);       // This skips i to almost the end of the tree, any remaining entries
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

   return 0;
}


void IncSpectra()
{                               // For testing only
   //hTestSpectrum->SetBinContent(1,1000);
}



// Function to parse Mnemonic name:
void ParseMnemonic(std::string * name, Mnemonic * mnemonic)
{
   std::string buf;
   mnemonic->system.assign(*name, 0, 2);
   mnemonic->subsystem.assign(*name, 2, 1);
   buf.clear();
   buf.assign(*name, 3, 2);
   mnemonic->arrayposition = atoi(buf.c_str());
   mnemonic->arraysubposition.assign(*name, 5, 1);
   mnemonic->collectedcharge.assign(*name, 6, 1);
   buf.clear();
   buf.assign(*name, 7, 2);
   mnemonic->segment = atoi(buf.c_str());
   mnemonic->outputsensor.assign(*name, 9, 1);
}

int Col2Num(char Colour)
{
   switch (Colour) {
   case 'B':
      return 0;
   case 'G':
      return 1;
   case 'R':
      return 2;
   case 'W':
      return 3;
   default:
      return -1;
   }
}

char Num2Col(int Crystal)
{
   switch (Crystal) {
   case 0:
      return 'B';
   case 1:
      return 'G';
   case 2:
      return 'R';
   case 3:
      return 'W';
   default:
      return 'X';
   }
}

int ReadCalibrationFile(std::string filename, vector < string > *CalibNames, vector < vector < float >>*CalibValues)
{

   if (Config.PrintBasic) {
      printf("Reading calibration file %s...\t", filename.c_str());
   }
   ifstream file;
   file.open(filename);
   if (!file) {
      if (Config.PrintBasic) {
         printf("could not open file.\n");
      }
      return -1;
   } else {
      if (Config.PrintBasic) {
         printf("File opened.\n");
      }
   }

   std::string line;
   char name[CHAR_BUFFER_SIZE], temp[CHAR_BUFFER_SIZE];
   float g0, g1, g2;
   int n = 0;

   while (getline(file, line)) {        // loop lines in file
      if (line.c_str()[0] != '#') {     // skip commented lines
         g0 = 0.0;              // reset values
         g1 = 0.0;
         g2 = 0.0;
         sscanf(line.c_str(), "%s %s %f %f %f", name, temp, &g0, &g1, &g2);     // grab name and 3 gains
         vector < float >GainTemp;
         GainTemp.push_back(g0);
         GainTemp.push_back(g1);
         GainTemp.push_back(g2);
         CalibValues->push_back(GainTemp);      // store values
         CalibNames->push_back(name);   // store name
         n += 1;                // count
      }
   }
   return n;
}


float CalibrateEnergy(int Charge, std::vector < float >Coefficients)
{

   float ChargeF = (float) Charge + rand1.Uniform();
   float TempInt = 125.0;
   float Energy = 0.0;
   if (Coefficients.size() == 0) {
      return ChargeF;
   }
   if (INTEGRATION != 0) {
      TempInt = INTEGRATION;
   }
   for (int i = 0; i < Coefficients.size(); i++) {
      Energy += Coefficients[i] * pow((ChargeF / TempInt), i);
   }
   return Energy;
}


float CalibrateWaveEnergy(float Charge, std::vector < float >Coefficients)
{

   Charge = Charge + rand1.Uniform();
   float Energy = 0.0;
   if (Coefficients.size() == 0) {
      return Charge;
   }
   for (int i = 0; i < Coefficients.size(); i++) {
      Energy += Coefficients[i] * pow(Charge, i);
   }
   return Energy;
}



float CalcWaveCharge(std::vector < int >wavebuffer)
{
   int Samp, Length;
   float Charge = 0.0;
   float Initial = 0.0;
   float Final = 0.0;
   Length = wavebuffer.size();

   if (wavebuffer.size() < Config.WaveInitialSamples + Config.WaveFinalSamples) {
      return 0.0;               // return 0 if wave too short
   }
   for (Samp = 0; Samp < Config.WaveInitialSamples; Samp++) {
      Initial += wavebuffer.at(Samp);
   }
   Initial /= Config.WaveInitialSamples;
   for (Samp = 0; Samp < Config.WaveFinalSamples; Samp++) {
      Final += wavebuffer.at(Length - Samp - 1);        // -1 because final sample in wbuffer seems to be spurious
   }
   Final /= Config.WaveFinalSamples;
   Charge = Final - Initial;
   return Charge;
}







