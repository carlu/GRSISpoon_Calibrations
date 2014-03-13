//To compile:
// g++ Main.C CoincEff.C Calib.C PropXtalk.C CalibTools.C CalibOffline.C -I$GRSISYS/include --std=c++0x -o Sort $GRSISYS/libraries/TigFormat/libFormat.so $GRSISYS/libraries/libCalManager.so $GRSISYS/libraries/libRootIOManager.so -O2 `root-config --cflags --libs` -lTreePlayer -lSpectrum -lgsl -lgslcblas -g
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

// canvases for calibration.  To be shared between Calib() and CalibOffline() functions
TCanvas *cCalib1, *cCalib1a, *cCalib2, *cWave1, *ctemp;

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
void PrintHelp();

void SortTree(const char *fn);
void IncSpectra();
int ReadCalibrationFile(std::string filename, vector < string > *EnCalibNames,vector < vector < float >> *EnCalibValues);

void CoincEff(std::vector < TTigFragment > &ev);
void InitCoincEff();
void FinalCoincEff();

void Calib(std::vector < TTigFragment > &ev);
void InitCalib();
void FinalCalib();

int CalibOffline(std::string filename);

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
   
   // Offline calibration here
   // need to initialise TCanvas's for Calib and CalibOffline here:
   if(Config.RunOffCal == 1 || Config.RunCalibration == 1) {
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
   if(Config.RunOffCal== 1) {
      for(i=0;i<Config.files.size();i++) {
         if(strncmp(Config.files.at(i).c_str(),Config.CalOut.c_str(),8)==0) {
            cout << "Attempting offline calibration on histograms in file: " << Config.files.at(i) << endl;
            CalibOffline(Config.files.at(i));
         }
      }
   
      return 0;
   }
   // ---- End of offline calibration -----

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

   for (ChainEvent = 0; ChainEvent < NumChainEntries; ChainEvent++) {

      //cout << "HERE!" << endl;

      Chain->LoadTree(ChainEvent);
      TreeNum = Chain->GetTreeNumber();

      if (TreeNum != LastTreeNum) {
         cout << "Switching to TreeNum " << TreeNum + 1 << " from " << LastTreeNum +
             1 << " at chain entry " << ChainEvent << endl;
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

      for (int TreeEvent = 0; TreeEvent < NumTreeEvents; TreeEvent++) {
         evFrags.clear();
         int FragNum = 1;

         while (Tree->GetEntryWithIndex(TreeEvent, FragNum++) != -1) {
            evFrags.push_back(*pFrag);
            FragCount++;
            if (DEBUG_TREE_LOOP) {
               cout << "FragNum: " << FragNum << " j: " << TreeEvent;
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

         // do something with the evFrag vector (contains a built event)... ProcessEvent(evFrags); 
         if (Config.EventLimit > 0 && EventCount >= Config.EventLimit) {
            cout << "Maximum number of events (" << Config.EventLimit << ") reached.  Terminating..." << endl;
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
         if (PRINT_OUTPUT && (TreeEvent % PRINT_FREQ) == 0) {
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

      if (Config.EventLimit > 0 && EventCount >= Config.EventLimit) {
         cout << "Maximum number of events (" << Config.EventLimit << ") reached.  Terminating..." << endl;
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

   while (getline(file, line)) {       // loop lines in file
      if (line.c_str()[0] != '#') {     // skip commented lines
         g0 = 0.0;   // reset values
         g1 = 0.0;
         g2 = 0.0;
         sscanf(line.c_str(), "%s %f %f %f", name, &g0, &g1, &g2);  // grab name and 3 gains
         vector < float >GainTemp;  
         GainTemp.push_back(g0);
         GainTemp.push_back(g1);
         GainTemp.push_back(g2);
         CalibValues->push_back(GainTemp);   // store values
         CalibNames->push_back(name);        // store name
         n += 1;  // count
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
   Config.RunOffCal = SORT_OFFCAL;
   Config.RunEfficiency = SORT_EFF;
   Config.RunPropCrosstalk = SORT_PROP;
   Config.RunWaveform = SORT_WAVES;
   Config.RunDiffCrosstalk = SORT_DIFF;
   
   Config.OutPath  = "./";
   Config.CalOut   = "CalibOut.root";
   Config.CalOfOut = "CalibOffOut.root";
   Config.EffOut   = "CoincEffOut.root";
   Config.EffTxtOut= "CoincEffOut.txt";
   Config.PropOut  = "PropXTalkOut.root";
     
   Config.PrintBasic = PRINT_OUTPUT;
   Config.PrintFrequency = PRINT_FREQ;
   Config.PrintVerbose = PRINT_VERBOSE;
   
   Config.EventLimit = MAX_EVENTS;
   
   Config.EnergyCalibrationFile = "./ECal.txt";
   Config.HaveAltEnergyCalibration = 0;
   Config.WaveCalibrationFile = "./WCal.txt";
   Config.HaveWaveCalibration = 0;
   
   float Sources[3][10] = {
      {1173.237, 1332.501},  // 60Co
      {121.7817, 1408.006, 244.6975, 344.2785, 411.116, 778.9040, 964.079, 1112.074},  // 152Eu
      {344.2785, 1408.006, 244.6975, 411.116, 778.9040, 964.079, 1112.074}   // 152Eu (no 121)
   };
   
   // Global physics settings
   // ------------------------------------------
   // Source information
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
   // Properties of waveforms stored in the data
   Config.WaveformSamples = WAVE_SAMPS;
   Config.WaveInitialSamples = INITIAL_SAMPS;
   Config.WaveFinalSamples = FINAL_SAMPS;
   // Thresholds
   Config.EnergyThresh = 5; // keV
   Config.ChargeThresh = 100;
   
   // Optons for calibration
   // Plots
   Config.PlotFits = 0;
   Config.PlotCalib = 0;
   Config.PlotCalibSummary = 0;
   
   Config.CalEnergy = FIT_EN;
   Config.CalWave   = FIT_WAVE_EN;
   Config.CalReport = OUTPUT_REPORT;
   
   Config.WriteFits = WRITE_FITS;
   
   
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
   // --calof : run calibration on spectrum file rather than fragment tree 
   // --eff : run efficiency
   // --prop: run propxtalk
   
   int i, j,n;
   bool test;
   bool RunConfGiven=0;
   
   /*cout << "List of arguments (" << argc << " total)" << endl;
   for(i=0;i<argc;i++) {
      cout << i << ": " << argv[i] << "\t";
   }
   cout << endl << endl;*/
   
   if(argc < 3) {
      cout << "No input file provided!" << endl << endl;
      PrintHelp();
      return -1;
   }
   
   for(i = 0; i < argc; i++) { // loop all args 
      cout << "i=" << i << endl;
      // Print help information
      if(strncmp(argv[i],"-h",2)==0) {
         PrintHelp();
         return -1;
      }
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
         if(strncmp(argv[i+1],"60Co",4)==0 || strncmp(argv[i+1],"Co60",4)==0 ||
              strncmp(argv[i+1],"60co",4)==0 || strncmp(argv[i+1],"co60",4)==0 ) {  // is it 60Co
            Config.SourceNumCore = 0;
            Config.SourceNumFront = 0;
            Config.SourceNumBack = 0;
            test=1;
         }
         if(strncmp(argv[i+1],"152Eu",5)==0 || strncmp(argv[i+1],"Eu152",5)==0 ||
              strncmp(argv[i+1],"152eu",5)==0 || strncmp(argv[i+1],"eu152",5)==0) {  // or is it 152Eu
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
      
      // Maximum number of events to process
      // -------------------------------------------
      if(strncmp(argv[i],"-n",2)==0) { 
         if(i>=argc-1 || strncmp(argv[i+1],"-",1)==0) {  // return error if no source given
            cout << "No number specified after \"-n\" option. (maximum number of events)" << endl;
            return -1;  
         }
         Config.EventLimit = atoi(argv[++i]);
         if(Config.EventLimit < 0) {
            cout << "Negative number of events specified.  What does that even mean?" << endl;
            return -1;
         }
      }
      
      // Maximum number of events to process
      // -------------------------------------------
      if(strncmp(argv[i],"-o",2)==0) { 
         if(i>=argc-1 || strncmp(argv[i+1],"-",1)==0) {  // return error if no source given
            cout << "No path specified after \"-o\" option. (output path)" << endl;
            return -1;  
         }
         Config.OutPath = argv[++i];
      }
      // Run options
      // -------------------------------------------
      if(strncmp(argv[i],"--",2)==0) {// Long option used.
         if(RunConfGiven==0) { // If run configuration given, clear defaults
            Config.RunCalibration = 0;
            Config.RunEfficiency = 0;
            Config.RunPropCrosstalk = 0;
            Config.RunWaveform = 0;
            Config.RunDiffCrosstalk = 0;
            RunConfGiven = 1;
         }
         // offline calibration
         if(strncmp(argv[i],"--calof",7)==0) {
            Config.RunOffCal = 1;
         }
         // calibration
         if(strncmp(argv[i],"--cal",5)==0) {
            Config.RunCalibration = 1;
         }
         // efficiency
         if(strncmp(argv[i],"--eff",5)==0) {
            Config.RunEfficiency = 1;
         }
         // proportional crosstalk
         if(strncmp(argv[i],"--prop",6)==0) {
            Config.RunPropCrosstalk = 1;
         }
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

void PrintHelp() {
   cout << "You seem confused, perhaps this will help.." << endl << endl;
   cout << "To run:" << endl << "./Sort -f InFile1 [InFile2...] " << endl;
   cout << "Options: " << endl;
   cout << "[-e (energy Calibration File)] - select alternate energy calibration file.  In the absense of an entry in this file, all channels will default to using the calibrated energy from the input TTree." << endl;
   cout << "[-w (Wave calibration file)] - select calibration for energy derived from waveforms.  No defaults. Required for succesfully running XTalk analysis." << endl;
   cout << "[-s (Source e.g. 60Co, 152eu)] - select source to be used for calibration." << endl;
   cout << "[-n N] - limit the number of events processed to N" << endl;
   cout << "[-o (path)] - save all output to (path) rather than ./" << endl;
   cout << "[--cal/--eff/--prop] - run the calibration, efficiency, or proportianal crosstalk parts of the code on a fragment tree input." << endl;
   cout << "[--calof] - run offline calibration on a histogram file (CalibOutXXXX.root).  This option overrides all other run options." << endl; 
   cout << endl;
}






