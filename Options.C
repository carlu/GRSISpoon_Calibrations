// Code related to run configuration and command line options
// -------------------------------------------------------------

// C/C++ libraries:
using namespace std;
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string.h>

#include "Options.h"

// Functions
int LoadDefaultSettings();
int ReadCommandLineSettings(int argc, char **argv);
int ReadConfigFile(std::string filename);
void PrintHelp();

// Storing run settings
RunConfig Config;

// Load sensible defaults before reading any options
int LoadDefaultSettings()
{

   Config.RunCalibration = 0;
   Config.RunSpecCal = 0;
   Config.RunEfficiency = 0;
   Config.RunSpecEff = 0;
   Config.RunPropCrosstalk = 0;
   Config.RunWaveform = 0;
   Config.RunDiffCrosstalk = 0;

   Config.OutPath = "./";

   Config.PrintBasic = 1;
   Config.PrintFrequency = PRINT_FREQ;
   Config.PrintVerbose = 0;

   Config.EventLimit = MAX_EVENTS;

   Config.ROOT_MaxVirtSize = ROOT_VIRT_SIZE;

   Config.EnergyCalibrationFile = "./ECal.txt";
   Config.HaveAltEnergyCalibration = 0;
   Config.WaveCalibrationFile = "./WCal.txt";
   Config.HaveWaveCalibration = 0;
   
   Config.ConfigFile = "./Config.txt";

   // Global physics settings
   // ------------------------------------------

   // Source information
   float Sources[4][10] = {
      {1173.237, 1332.501},     // 60Co
      {121.7817, 1408.006, 244.6975, 344.2785, 411.116, 778.9040, 964.079, 1112.074},   // 152Eu
      {344.2785, 1408.006, 244.6975, 411.116, 778.9040, 964.079, 1112.074},     // 152Eu (no 121)
      {276.398, 356.017, 80.9971, 302.853, 383.851}     // 133Ba
   };
   // Push source lines to vector of vectors.
   // This is horrible.  Need to find another way to initialise these vectors.
   vector < float >SourceTemp;
   for (int i = 0; i < 2; i++) {
      SourceTemp.push_back(Sources[0][i]);
   }
   Config.Sources.push_back(SourceTemp);
   SourceTemp.clear();
   for (int i = 0; i < 7; i++) {
      SourceTemp.push_back(Sources[1][i]);
   }
   Config.Sources.push_back(SourceTemp);
   SourceTemp.clear();
   for (int i = 0; i < 7; i++) {
      SourceTemp.push_back(Sources[2][i]);
   }
   Config.Sources.push_back(SourceTemp);
   SourceTemp.clear();
   for (int i = 0; i < 5; i++) {
      SourceTemp.push_back(Sources[3][i]);
   }
   Config.Sources.push_back(SourceTemp);

   // Properties of waveforms stored in the data
   Config.WaveformSamples = 200;
   Config.WaveInitialSamples = 65;
   Config.WaveFinalSamples = 65;
   // Thresholds
   Config.EnergyThresh = 5;     // keV
   Config.ChargeThresh = 100;
   // Array mode
   Config.HighEffMode = 0;

   // Options for Calib() and CalibOffline()
   // What to do
   Config.CalEnergy = 1;        // Calibrate charge spectra?
   Config.CalWave = 1;          // Calibrate wave spectra?
   Config.CalReport = 1;        // Generate report on fits etc
   Config.CalFile = 0;          // Generate .cal file as used by GRSISpoon 
   memset(&Config.CalList, 1, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
   Config.CalListProvided = 0;
   // Output
   Config.CalOut = "CalibOut.root";     // Name for Tree calibration output file
   Config.CalName = Config.CalOut.substr(0, Config.CalOut.size() - 5);  // name expected for hist calib input.  
   // .root stripped because run number may exist 
   Config.AnaName = "his";
   Config.CalSpecOut = "CalibSpecOut.root";
   Config.WriteFits = 1;
   // Plots
   Config.PlotFits = 0;
   Config.PlotCalib = 0;
   Config.PlotResidual = 0;
   Config.PlotCalibSummary = 0;
   memset(&Config.CalibPlots, 0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
   // Calibrating options
   Config.FitZero = 0;
   memset(&Config.ManualPeakSelect, 0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
   Config.ManualPeakCorrection = 1;

   // Options for CoincEff()
   Config.EffOut = "CoincEffOut.root";
   Config.EffTxtOut = "CoincEffOut.txt";

   Config.EffSimRefFileName = "SimEffRef.txt";
   Config.EffExpRefFileName = "ExpEffRef.txt";

   Config.OutputEff = 1;
   // Plots
   Config.PlotEff = 1;

   // Options for PropXtakl()
   Config.PropOut = "PropXTalkOut.root";

   return 0;
}

// Override defaults with options from commandd line. 
// Also some options here are required i.e. file and function to perform
int ReadCommandLineSettings(int argc, char **argv)
{

   // -f : input (f)iles
   // -e : (E)nergyCalFile
   // -w : (W)aveCalFile
   // -s : specify (s)ource 60Co 152Eu...
   // -c : specify (c)onfig options file  #COMMENT\nNAME VALUE\nNAME VALUE
   // -o : (o)utput path (prepended to all output files)
   // -h : print (h)elp and exit
   // -vr: Set ROOT virtual memory size in bytes
   // -v : (v)erbose
   // -q : (Q)uiet
   // -n : max (n)umber of events

   // -p : (p)lot (Clover) (Crystal) (Seg)
   // -mp: (m)anual (p)eak  (Clover) (Crystal) (Seg)  : manually selcect peaks to be used on this seg
   //                                                    or all segs if none specified

   // -z : add extra calibration point at (z)ero  i.e. 0ch = 0keV
   // -d : select (d)etector to be calibrated.  

   // --cal : run calibration
   // --calspec : run calibration on spectrum file rather than fragment tree 
   // --eff : run efficiency
   // --prop: run propxtalk

   int i, n;
   unsigned int FileNum;
   int Plot[3];
   int FitList[3];
   int Limits[3] = { CLOVERS, CRYSTALS, SEGS + 2 };
   bool test;
   bool RunConfGiven = 0;       // used to record if required operation has been specified.  
   bool ConfigFileLoaded = 0;

   if (argc < 3) {
      cout << "No input file provided!" << endl << endl;
      PrintHelp();
      return -1;
   }

   for (i = 0; i < argc; i++) { // loop all args 
      // Print help information
      if (strncmp(argv[i], "-h", 2) == 0) {
         PrintHelp();
         return -1;
      }
      // Run files
      // -------------------------------------------
      if (strncmp(argv[i], "-f", 2) == 0) {     // if option is input file
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no file
            cout << "No file specified after \"-f\" option." << endl;
            return -1;
         }
         while (strncmp(argv[i + 1], "-", 1) > 0) {     // loop files 
            Config.files.push_back(argv[++i]);  // add to vector and increment i
            if (i >= argc - 1) {
               break;           // break if at last item in arg list
            }
         }
         if (Config.PrintBasic) {
            cout << "Input files:  " << endl;
            for (FileNum = 0; FileNum < Config.files.size(); FileNum++) {       // print list of files back to screen
               cout << "\t" << Config.files.at(FileNum) << endl;
            }
         }
      }
      // Load Configuration file
      // -------------------------------------------
      if (strncmp(argv[i], "-c", 2) == 0) {     // if option is Ecal file
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no file
            cout << "No file specified after \"-c\" option." << endl;
            return -1;
         }
         if(ReadConfigFile(argv[++i]) == 0) {
            ConfigFileLoaded = 1;
            Config.ConfigFile = argv[i];
            cout << "Configuration loaded from " << Config.ConfigFile << endl;
         }
         else {
            cout << "Failed to load configuration from " << Config.ConfigFile << endl;
         }
         
      }
      // Energy Calibration file
      // -------------------------------------------
      if (strncmp(argv[i], "-e", 2) == 0) {     // if option is Ecal file
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no file
            cout << "No file specified after \"-e\" option." << endl;
            return -1;
         }
         Config.EnergyCalibrationFile = argv[++i];      // otherwise set filename
         Config.HaveAltEnergyCalibration = 1;
         cout << "Energy calibration to be loaded from " << Config.EnergyCalibrationFile << endl;
      }
      // Waveform Calibration file
      // -------------------------------------------
      if (strncmp(argv[i], "-w", 2) == 0) {     // if option is Ecal file
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no file
            cout << "No file specified after \"-w\" option." << endl;
            return -1;
         }
         Config.WaveCalibrationFile = argv[++i];        // otherwise set filename
         Config.HaveWaveCalibration = 1;
         cout << "Wave calibration to be loaded from " << Config.WaveCalibrationFile << endl;
      }
      // Source specification
      // -------------------------------------------
      if (strncmp(argv[i], "-s", 2) == 0) {     // if option is source spec
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no source given
            cout << "No source specified after \"-s\" option." << endl;
            return -1;
         }
         while (strncmp(argv[i + 1], "-", 1) > 0) {     // loop files 
            test = 0;
            if (strncmp(argv[i + 1], "60Co", 4) == 0 || strncmp(argv[i + 1], "Co60", 4) == 0 || strncmp(argv[i + 1], "60co", 4) == 0 || strncmp(argv[i + 1], "co60", 4) == 0) { // is it 60Co
               Config.SourceNumCore.push_back(0);
               Config.SourceNumFront.push_back(0);
               Config.SourceNumBack.push_back(0);
               test = 1;
            }
            if (strncmp(argv[i + 1], "152Eu", 5) == 0 || strncmp(argv[i + 1], "Eu152", 5) == 0 || strncmp(argv[i + 1], "152eu", 5) == 0 || strncmp(argv[i + 1], "eu152", 5) == 0) {     // or is it 152Eu
               Config.SourceNumCore.push_back(1);
               Config.SourceNumFront.push_back(1);
               Config.SourceNumBack.push_back(2);
               test = 1;
            }
            if (strncmp(argv[i + 1], "133Ba", 5) == 0 || strncmp(argv[i + 1], "Ba133", 5) == 0 || strncmp(argv[i + 1], "133ba", 5) == 0 || strncmp(argv[i + 1], "ba133", 5) == 0) {     // or is it 133Ba
               Config.SourceNumCore.push_back(3);
               Config.SourceNumFront.push_back(3);
               Config.SourceNumBack.push_back(3);
               test = 1;
            }
            if (test == 0) {    // or is it somethimng else
               cout << "Source not recognised: " << argv[i + 1] << endl;
               return -1;
            }
            i += 1;
            if (i >= argc - 1) {
               break;           // break if at last item in arg list
            }
         }




      }
      // Maximum number of events to process
      // -------------------------------------------
      if (strncmp(argv[i], "-n", 2) == 0) {
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no source given
            cout << "No number specified after \"-n\" option. (maximum number of events)" << endl;
            return -1;
         }
         Config.EventLimit = atoi(argv[++i]);
         if (Config.EventLimit < 0) {
            cout << "Negative number of events specified.  What does that even mean?" << endl;
            return -1;
         }
         if (Config.EventLimit == 0) {
            cout << "No event limit set." << endl;
         } else {
            cout << "Run limited to " << Config.EventLimit << " events." << endl;
         }
      }
      // Maximum number of events to process
      // -------------------------------------------
      if (strncmp(argv[i], "-o", 2) == 0) {
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no source given
            cout << "No path specified after \"-o\" option. (output path)" << endl;
            return -1;
         }
         Config.OutPath = argv[++i];
      }
      // ROOT virtual memory max size
      // -----------------------------------
      if (strncmp(argv[i], "-vr", 3) == 0) {
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {
            cout << "No size specified after \"-vr\" option" << endl;
            return -1;
         }
         if (Config.EventLimit < 0) {
            cout << "Negative Virtual RAM size specified" << endl;
            return -1;
         } else {
            Config.ROOT_MaxVirtSize = atoi(argv[++i]) * 1024u * 1024u;
            cout << "ROOT TTree->SetMaxVirtualSize( " << Config.ROOT_MaxVirtSize << " )" << endl;
         }

      }
      // Verbose mode
      // -------------------------------------------
      if (strncmp(argv[i], "-v", 2) == 0) {
         Config.PrintBasic = 1;
         Config.PrintVerbose = 1;
      }
      // Quiet mode   (clearly verbose and quiet will override eachother)
      // -------------------------------------------
      if (strncmp(argv[i], "-q", 2) == 0) {
         Config.PrintBasic = 0;
         Config.PrintVerbose = 0;
      }
      // Plots 
      // ----------------------------------------
      if (strncmp(argv[i], "-p", 2) == 0) {
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no file
            cout << "Plotting all fits." << endl;
            Config.PlotFits = 1;
            memset(&Config.CalibPlots, 1, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
         } else {
            Config.PlotFits = 1;
            n = 0;
            while (strncmp(argv[i + 1], "-", 1) > 0) {  // loop plot items 
               Plot[n] = atoi(argv[++i]);
               cout << "n,Plot[n] = " << n << ", " << Plot[n] << endl;
               if ((Plot[n] < 0) || (Plot[n] > Limits[n])) {
                  cout << "Error with plot specification!" << endl;
                  return -1;
               }
               n++;
               if (n == 3) {
                  break;
               }
               if (i >= argc - 1) {
                  cout << "Need to specify (Clover) (Crystal) (Seg) for plots" << endl;
                  return -1;
               }
            }
            cout << "Plotting Cl: " << Plot[0] << " Cr: " << Plot[1] << " Seg: " << Plot[2] << endl;
            Config.CalibPlots[Plot[0] - 1][Plot[1]][Plot[2]] = 1;
         }
      }
      // Manual PEak Selection
      // -------------------------------------------
      if (strncmp(argv[i], "-mp", 3) == 0) {
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // if last arg or next is new option
            cout << "--------------------------------------" << endl;
            cout << "Manual Peak Selection For Calibration!" << endl;
            cout << "--------------------------------------" << endl;
            memset(&Config.ManualPeakSelect, 1, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
         } else {
            n = 0;
            while (strncmp(argv[i + 1], "-", 1) > 0) {  // loop plot items 
               Plot[n] = atoi(argv[++i]);
               cout << "n,Plot[n] = " << n << ", " << Plot[n] << endl;
               if ((Plot[n] < 0) || (Plot[n] > Limits[n])) {    // allowed range is the same as for plotting
                  cout << "Error with specification of manual peak select!" << endl;
                  return -1;
               }
               n++;
               if (n == 3) {
                  break;
               }
               if (i >= argc - 1) {
                  cout << "Need to specify (Clover) (Crystal) (Seg) for manual peak select" << endl;
                  return -1;
               }
            }
            cout << "Manually selecting peaks for Cl: " << Plot[0] << " Cr: " << Plot[1] << " Seg: " << Plot[2] << endl;
            Config.ManualPeakSelect[Plot[0] - 1][Plot[1]][Plot[2]] = 1;
            Config.CalibPlots[Plot[0] - 1][Plot[1]][Plot[2]] = 1;       // Also set this channel to plot so we can see it!
         }
      }
      // add 0ch = 0keV to calibration
      // -------------------------------------------
      if (strncmp(argv[i], "-z", 2) == 0) {
         Config.FitZero = 1;
      }
      // List of detectors to fit
      // -------------------------------------------
      if (strncmp(argv[i], "-d", 2) == 0) {
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {
            cout << "Need to specify (Clover) (Crystal) (Seg) with \"-d\" option." << endl;
            return -1;
         } else {
            n = 0;
            while (strncmp(argv[i + 1], "-", 1) > 0) {
               FitList[n] = atoi(argv[++i]);
               if ((FitList[n] < 0) || (FitList[n] > Limits[n])) {
                  cout << "Error with plot specification!" << endl;
                  return -1;
               }
               n++;
               if (n == 3) {
                  break;
               }
               if (i >= argc - 1) {
                  cout << "Need to specify (Clover) (Crystal) (Seg) with \"-d\" option." << endl;
                  return -1;
               }
            }
            if (Config.CalListProvided == 0) {
               Config.CalListProvided = 1;
               memset(&Config.CalList, 0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
            }
            Config.CalList[FitList[0] - 1][FitList[1]][FitList[2]] = 1;
         }
      }
      // Run options
      // -------------------------------------------
      if (strncmp(argv[i], "--", 2) == 0) {     // Long option used.
         if (RunConfGiven == 0) {       // If run configuration given, clear defaults
            Config.RunCalibration = 0;
            Config.RunEfficiency = 0;
            Config.RunPropCrosstalk = 0;
            Config.RunWaveform = 0;
            Config.RunDiffCrosstalk = 0;
            RunConfGiven = 1;
         }
         // offline calibration
         if (strncmp(argv[i], "--calspec", 7) == 0) {
            Config.RunSpecCal = 1;
         }
         // calibration
         if (strncmp(argv[i], "--cal", 5) == 0) {
            Config.RunCalibration = 1;
         }
         // efficiency
         if (strncmp(argv[i], "--eff", 5) == 0) {
            Config.RunEfficiency = 1;
         }
         // proportional crosstalk
         if (strncmp(argv[i], "--prop", 6) == 0) {
            Config.RunPropCrosstalk = 1;
         }
      }
   }
   // Look for default configuration file if non given
   if(ConfigFileLoaded == 0) {
      
   }   

   return 0;
}


int ReadConfigFile(std::string filename)        // this is a dummy function
{                               // One day it should read in a config text file

   cout << filename << endl;
   return 0;
}



void PrintHelp()
{
   cout << "You seem confused, perhaps this will help.." << endl << endl;
   cout << "To run:" << endl << "./Sort[Histos/Trees] -f InFile1 [InFile2...] " << endl;
   cout << "Options: " << endl << endl;
   cout <<
       "\t[-e (energy Calibration File)] - select alternate (e)nergy calibration file.  In the absense of an entry in this file, all channels will default to using the calibrated energy from the input TTree."
       << endl;
   cout <<
       "\tFormat for calibration file should be (TIGNOM) (some other string e.g. chg or wavechg) (g0) (g1) [(g2) (g3)].  Lines begining \"#\" will be ignored."
       << endl << endl;
   cout <<
       "[-w (Wave calibration file)] - select calibration for energy derived from (w)aveforms.  No defaults. Required for succesfully running XTalk analysis. Format same as for energy calibration file."
       << endl << endl;
   cout << "[-s (Source e.g. 60Co, 152eu, Ba133)] - select (s)ource to be used for calibration." << endl;
   cout << "\tThis is required to run calibration." << endl << endl;
   cout << "[-n N] - limit the (n)umber of events processed to N" << endl << endl;
   cout << "[-o (path)] - save all (o)utput to (path) rather than ./" << endl << endl;
   cout <<
       "[--cal/--eff/--prop] - run the calibration, efficiency, or proportianal crosstalk parts of the code on a fragment tree input."
       << endl << endl;
   cout <<
       "[--calspec] - run calibration on a histogram file (CalibOutXXXX.root, hisXXXX.root).  This option overrides all other run options."
       << endl;
   cout << "\tCalibOutX.root format assumes spectra named according to those produced by --cal." << endl;
   cout <<
       "\thisX.root assumes no wave spectra and format ChrgXXXX names with numbers according to TIGRESS DAQ analyser convention"
       << endl << endl;
   cout << "[-z] - Include an additional point at 0ch = 0keV in any calibrations." << endl << endl;
   cout <<
       "[-q] - Quiet mode.  [-v] - Verbose mode.  Default prints progress through run and configuration.  Quiet doesn't.  Verbose also prints other stuff"
       << endl << endl;
   cout <<
       "[-p (Clover) (Crystal) (Seg)] - Plot certain segment when doing calibration.  With no additional arguments will plot all."
       << endl << endl;

   cout <<
       "[-d (Clover) (Crystal) (Seg)] - Selects a particular (d)etection element to fit.  Defaults to all. Affect --calspec only"
       << endl << endl;
   cout <<
       "[-mp (Clover) (Crystal) (Seg)] - Manually selects calibration peaks for a particular (d)etection element to   Defaults to all."
       << endl << endl;
   cout << "[-vr (VRam in Mb)] - Specify in Mb what should be passed to Tree->SetMaxVirtualSize()." << endl;
   cout << "\tOnly effects SortTrees, not SortHistos." << endl;
   cout << "\tToo small and sort will be slow, too large and it will crash." << endl;
   cout << "\tIt is printed back in bytes.  (Overflows above 4095Mb, check if unsure)" << endl;

   cout << endl;
}
