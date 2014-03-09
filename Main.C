//To compile:
// g++ Main.C CoincEff.C Calib.C PropXtalk.C CalibTools.C -I$GRSISYS/include --std=c++0x -o Sort $GRSISYS/libraries/TigFormat/libFormat.so $GRSISYS/libraries/libCalManager.so $GRSISYS/libraries/libRootIOManager.so -O2 `root-config --cflags --libs` -lTreePlayer -lSpectrum -lgsl -lgslcblas -g
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

// TriScope libraries
#include "TTigFragment.h"
//#include "TFSPC_Info.h"
//#include "TSharc.h"
//#include "TTigress.h"
//#include "TRf.h"
//#include "TTriFoil.h"

// My libraries
// #include "CoincEff.h"
#include "Main.h"
//#include "Gains.h"

//#include "CalibrationManager.h"

// For tracking real time
TStopwatch watch;

TApplication *App;              // Pointer to root environment for plotting etc

// Rand for gain matching
static TRandom3 rand1;

// Counters
int FileCount = 0;
int EventCount = 0;
int FragCount = 0;
int BadEventCount = 0;

// Storing run settings
RunConfig Config;

// Storing alternate Calibration
vector < string > EnCalibNames;
vector < vector < float >> EnCalibValues;
// Waveform calibration
vector < string > WaveCalibNames;
vector < vector < float >> WaveCalibValues;

// Functions
int LoadDefaultSettings();
int ReadCommandLineSettings(int argc, char **argv);

void SortTree(const char *fn);
void IncSpectra();
int ReadCalibrationFile(std::string filename, vector < string > *EnCalibNames,vector < vector < float >> *EnCalibValues);

void CoincEff(std::vector < TTigFragment > &ev);
void InitCoincEff();
void FinalCoincEff();

void Calib(std::vector < TTigFragment > &ev);
void InitCalib();
void FinalCalib();

void PropXtalk(std::vector < TTigFragment > &ev);
void InitPropXtalk();
void FinalPropXtalk();

int main(int argc, char **argv)
{
   // Variables, Constants, etc
   int i, j;

   // Set default and read custom options
   LoadDefaultSettings();
   for(i=0;i<Config.Sources.size();i++) {
      cout << "Source " << i << endl;
      for(j=0;j<Config.Sources.at(i).size();j++) {
         cout << "En " << j << " = " << Config.Sources.at(i).at(j) << endl;
      } 
   }
   if(ReadCommandLineSettings(argc, argv) < 0) {
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
      cout << "Alternate energy calibratrion values: " << endl;
      for (i = 0; i < NumCal; i++) {
         cout << i;
         cout << ": " << EnCalibNames.at(i);
         cout << " " << EnCalibValues.at(i)[0] << " " << EnCalibValues[i][1] << " " << EnCalibValues[i][2] << endl;
      }
   }
   // Load gain coefficients for waveforms
   if(Config.HaveWaveCalibration) {
      int NumCal;
      NumCal = ReadCalibrationFile(Config.WaveCalibrationFile, &WaveCalibNames, &WaveCalibValues);
      cout << "Wave energy calibratrion values: " << endl;
      //cout << WaveCalibNames.size() << endl << WaveCalibValues.size() << endl;
      for (i = 0; i < NumCal; i++) {
         cout << i;
         cout << ": " << WaveCalibNames.at(i);
         cout << " " << WaveCalibValues.at(i)[0] << " " << WaveCalibValues[i][1] << " " << WaveCalibValues[i][2] << endl;
      }
   }
   // Initialise spectra   

   if (Config.RunEfficiency) {
      cout << "Initialising Efficiency Spectra..." << endl;
      InitCoincEff();
   }
   if (Config.RunCalibration) {
      cout << "Initialising Calibration Spectra..." << endl;
      InitCalib();
   }
   if (Config.RunPropCrosstalk) {
      cout << "Initialising Cross-talk Spectra..." << endl;
      InitPropXtalk();
   }

   TChain *Chain = new TChain("FragmentTree");

   // Setup input file list                                        
   for (i = 0; i < Config.files.size(); i++) {
      Chain->Add(Config.files.at(i).c_str());
      FileCount++;
   }

   TTigFragment *pFrag = 0;
   std::vector < TTigFragment > evFrags;

   int nTrees = Chain->GetNtrees();
   int NumChainEntries = Chain->GetEntries();
   int NumChainEvents = (int) Chain->GetMaximum("TriggerId");   // This doesn't work, TrigID reset for each tree on chain.

   cout << "Chain Entries (frags) : " << NumChainEntries << endl;
   cout << "Chain Events          : " << NumChainEvents << endl;

   int TreeNum = -1;
   int LastTreeNum = -1;
   int NumTreeEntries = -1;
   int NumTreeEvents = -1;

   for (i = 0; i < NumChainEntries; i++) {

      //cout << "HERE!" << endl;

      Chain->LoadTree(i);
      TreeNum = Chain->GetTreeNumber();

      if (TreeNum != LastTreeNum) {
         cout << "Switching to TreeNum " << TreeNum + 1 << " from " << LastTreeNum +
             1 << " at chain entry " << i << endl;
         LastTreeNum = TreeNum;
      } else {
         continue;              // This works in conjuncion with the line below "i += (NumTreeEntries - 10);"
         // if i still falls in the same tree, then it is incremented without sorting until the next tree
      }

      //cout << "HERE!!" << endl;

      TTree *Tree = Chain->GetTree();

      if (!Tree->GetTreeIndex()) {
         printf("\nTree Index not found, Building index...");
         fflush(stdout);
         Tree->BuildIndex("TriggerId", "FragmentId");
         printf("  Done\n");
         fflush(stdout);
      }

      TBranch *Branch = Tree->GetBranch("TTigFragment");

      //cout << "HERE!!!" << endl;

      Branch->SetAddress(&pFrag);
      Tree->SetMaxVirtualSize(ROOT_VIRT_SIZE);
      Branch->LoadBaskets();

      NumTreeEntries = Tree->GetEntries();
      NumTreeEvents = (int) Tree->GetMaximum("TriggerId");

      // cout << "HERE!!!! " << NumChainEntries << " " << NumTreeEntries << " " << NumTreeEvents << endl;

      for (int j = 0; j < NumTreeEvents; j++) {
         evFrags.clear();
         int FragNum = 1;

         //cout << endl << "------------------------------" << endl;
         //cout << "FragNum: " << FragNum << " j: " << j << endl;


         while (Tree->GetEntryWithIndex(j, FragNum++) != -1) {
            evFrags.push_back(*pFrag);
            FragCount++;
            if (DEBUG_TREE_LOOP) {
               cout << "FragNum: " << FragNum << " j: " << j;
            }
         }
         if (DEBUG_TREE_LOOP) {
            cout << endl;
         }
         //cout << "HERE!!!!! " << endl;

         if (FragNum < 3) {     // Init to 1, incremented on first GetEntryWithIndex, so should be at least 3 if any frags found      
            BadEventCount++;
            //cout << "HERE!!!!! " << BadEventCount << endl;
            continue;
         } else {
            EventCount++;
         }

         //cout << "HERE!!!!!!!!!!" << endl;
         // do something with the evFrag vector (contains a built event)... ProcessEvent(evFrags); 
         if (MAX_EVENTS > 0 && EventCount >= MAX_EVENTS) {
            break;
         }
         if (SORT_EFF) {
            CoincEff(evFrags);
         }                      //passing vector of built events.
         if (SORT_CALIB) {
            Calib(evFrags);
         }
         if (SORT_PROP) {
            PropXtalk(evFrags);
         }
         if (DEBUG_TREE_LOOP) {
            cout << "ev Num =  " << j << ", ev.size() = " << evFrags.size() << endl;
         }
         // Print info to stdout
         if (PRINT_OUTPUT && (j % PRINT_FREQ) == 0) {
            cout << "----------------------------------------------------------" << endl;
            cout << "  Tree  " << TreeNum + 1 << " / " << nTrees << endl;
            cout << "  Event " << EventCount +
                BadEventCount << " / " << NumChainEvents << "   (" << EventCount << " good, " << BadEventCount <<
                " empty)" << endl;
            cout << "  Frag  " << FragCount << " / " << NumChainEntries << endl;
            cout << "\t  in " << StopWatch.RealTime() << " seconds." << endl;
            StopWatch.Continue();
            cout << "----------------------------------------------------------" << endl;
         }

      }

      if (MAX_EVENTS > 0 && EventCount >= MAX_EVENTS) {
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

int ReadCalibrationFile(std::string filename,vector < string > *CalibNames,vector < vector < float >> *CalibValues)
{

   printf("Reading calibration file %s...\t", filename.c_str());
   ifstream file;
   file.open(filename);
   if (!file) {
      printf("could not open file.\n");
      return -1;
   } else {
      printf("File opened.\n");
   }

   std::string line;
   char name[16];
   float g0, g1, g2;
   int n = 0;

   while (getline(file, line)) {
      //cout << line << endl;    
      if (line.c_str()[0] != '#') {     // skip commented lines
         g0 = 0.0;
         g1 = 0.0;
         g2 = 0.0;
         sscanf(line.c_str(), "%s %f %f %f", name, &g0, &g1, &g2);
         //cout << name << " " << g0 << " " << g1 << " " << g2 << endl; 
         vector < float >GainTemp;
         GainTemp.push_back(g0);
         GainTemp.push_back(g1);
         GainTemp.push_back(g2);
         CalibValues->push_back(GainTemp);
         CalibNames->push_back(name);
         n += 1;
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

   // Print Wave for debugging
   //cout << "Printing Wave (" << Length << " samples): " << endl;
   //for(Samp=0;Samp<Length;Samp++) {
   //   cout << wavebuffer.at(Samp) << " ";
   //}
   //cout << endl << endl;
   //cout << "Initial: " << Initial << endl;
   
   if(wavebuffer.size() < INITIAL_SAMPS + FINAL_SAMPS) {
      return 0.0;  // return 0 if wave too short
   }

   for (Samp = 0; Samp < INITIAL_SAMPS; Samp++) {
      Initial += wavebuffer.at(Samp);
      //cout << Initial << " ";
   }
   Initial /= INITIAL_SAMPS;
   //cout << "Initial: " << Initial;
   for (Samp = 0; Samp < FINAL_SAMPS; Samp++) {
      Final += wavebuffer.at(Length - Samp - 1);        // -1 because final sample in wbuffer seems to be spurious
   }
   Final /= FINAL_SAMPS;
   //cout << " Final: " << Final;
   Charge = Final - Initial;
   //cout << " Charge: " << Charge << endl;
   return Charge;
}


int LoadDefaultSettings() {

   Config.RunCalibration = SORT_CALIB;
   Config.RunEfficiency = SORT_EFF;
   Config.RunPropCrosstalk = SORT_PROP;
   Config.RunWaveform = SORT_WAVES;
   Config.RunDiffCrosstalk = SORT_DIFF;
   
   Config.PrintBasic = PRINT_OUTPUT;
   Config.PrintFrequency = PRINT_FREQ;
   Config.PrintVerbose = PRINT_VERBOSE;
   
   Config.EventLimit = MAX_EVENTS;

   Config.WaveformSamples = WAVE_SAMPS;
   Config.WaveInitialSamples = INITIAL_SAMPS;
   Config.WaveFinalSamples = FINAL_SAMPS;

   
   Config.EnergyCalibrationFile = "./ECal.txt";
   Config.HaveAltEnergyCalibration = 0;
   Config.WaveCalibrationFile = "./WCal.txt";
   Config.HaveWaveCalibration = 0;
   
   float Sources[3][10] = {
      {1173.237, 1332.501},
      {121.7817, 1408.006, 244.6975, 344.2785, 411.116, 778.9040, 964.079, 1112.074},
      {344.2785, 1408.006, 121.7817, 244.6975, 411.116, 778.9040, 964.079, 1112.074}
   };
   vector <float> SourceTemp;
   
   for(int i =0; i<2; i++) {
      SourceTemp.push_back(Sources[0][i]);
   }
   Config.Sources.push_back(SourceTemp);
   SourceTemp.clear();
   for(int i =0; i<7; i++) {
      SourceTemp.push_back(Sources[1][i]);
   }
   Config.Sources.push_back(SourceTemp);
   SourceTemp.clear();
   for(int i =0; i<7; i++) {
      SourceTemp.push_back(Sources[2][i]);
   }
   Config.Sources.push_back(SourceTemp);
   
   Config.SourceNumCore = SOURCE_NUM_CORE;
   Config.SourceNumFront = SOURCE_NUM_FRONT;
   Config.SourceNumBack = SOURCE_NUM_BACK;
   
   
   return 0;
   
   
   
}

int ReadCommandLineSettings(int argc, char **argv) {
   
   // -f : input files
   // -e : EnergyCalFile
   // -w : WaveCalFile
   // -s : specify source 60Co 152Eu...
   // -c : specify config options file  #COMMENT\nNAME VALUE\nNAME VALUE
   // -o : output path (prepended to all output files)
   // -h : print help and exit
   // -v : verbose
   // -n : max number of events
   
   // --cal : run calibration
   // --eff : run efficiency
   // --prop: run propxtalk
   
   int i, j,n;
   bool test;
   
   cout << "List of arguments (" << argc << " total)" << endl;
   for(i=0;i<argc;i++) {
      cout << i << ": " << argv[i] << "\t";
   }
   cout << endl << endl;
   
   if(argc < 3) {
      cout << "I at least need to know which file to process!" << endl;
      return -1;
   }
   
   for(i = 0; i < argc; i++) { // loop all args 
      cout << "i=" << i << endl;
      // Run files
      // -------------------------------------------
      if(strncmp(argv[i],"-f",2)==0) { // if option is input file
         if(i>=argc-1 || strncmp(argv[i+1],"-",1)==0) {  // return error if no file
            cout << "No file specified after \"-f\" option." << endl;
            return -1;  
         }
         while(strncmp(argv[i+1],"-",1)>0) {  // loop files 
            Config.files.push_back(argv[++i]);  // add to vector and increment i
            if(i>=argc-1) {  
               break;  // break if at last item in arg list
            }
         }
         cout << "Input files:" << endl;  
         for (j=0; j<Config.files.size(); j++) {  // print list of files back to screen
            cout << Config.files.at(j) << endl;
         }          
      }
      
      // Energy Calibration file
      // -------------------------------------------
      if(strncmp(argv[i],"-e",2)==0) { // if option is Ecal file
         if(i>=argc-1 || strncmp(argv[i+1],"-",1)==0) {  // return error if no file
            cout << "No file specified after \"-e\" option." << endl;
            return -1;  
         }
         Config.EnergyCalibrationFile = argv[++i]; // otherwise set filename
         Config.HaveAltEnergyCalibration = 1;
         cout << "Energy calibration to be loaded from " << Config.EnergyCalibrationFile << endl;
      }
      
      // Waveform Calibration file
      // -------------------------------------------
      if(strncmp(argv[i],"-w",2)==0) { // if option is Ecal file
         if(i>=argc-1 || strncmp(argv[i+1],"-",1)==0) {  // return error if no file
            cout << "No file specified after \"-w\" option." << endl;
            return -1;  
         }
         Config.WaveCalibrationFile = argv[++i]; // otherwise set filename
         Config.HaveWaveCalibration = 1;
         cout << "Wave calibration to be loaded from " << Config.WaveCalibrationFile << endl;
      }
      
      // Source specification
      // -------------------------------------------
      if(strncmp(argv[i],"-s",2)==0) { // if option is source spec
         if(i>=argc-1 || strncmp(argv[i+1],"-",1)==0) {  // return error if no source given
            cout << "No source specified after \"-s\" option." << endl;
            return -1;  
         }
         test=0;
         if(strncmp(argv[i+1],"60Co",4)==0 || strncmp(argv[i+1],"Co60",4)==0) {  // is it 60Co
            Config.SourceNumCore = 0;
            Config.SourceNumFront = 0;
            Config.SourceNumBack = 0;
            test=1;
         }
         if(strncmp(argv[i+1],"152Eu",5)==0 || strncmp(argv[i+1],"Eu152",5)==0) {  // or is it 152Eu
            Config.SourceNumCore = 1;
            Config.SourceNumFront = 1;
            Config.SourceNumBack = 2;
            test=1;
         }
         if(test==0) {  // or is it somethimng else
            cout << "Source not recognised!" << endl;
            return -1;
         }
         i += 1;
      }
      
      // General Configuration file
      // -------------------------------------------
      

   }
   
   
   return 0;
}

int ReadOptionsFile(std::string filename ) {
   cout << filename << endl;
   return 0;
}








