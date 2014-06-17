// Code related to run configuration and command line options
// -------------------------------------------------------------

// C/C++ libraries:
using namespace std;
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string.h>
#include <map>

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
   // Settings related to how the code works
   // ---------------------------------------
   // Which parts of code to run
   Config.RunCalibration = 0;
   Config.RunSpecCal = 0;
   Config.RunSegCoreCorrelation = 0;
   Config.RunEfficiency = 0;
   Config.RunSpecEff = 0;
   Config.RunPropCrosstalk = 0;
   Config.RunWaveform = 0;
   Config.RunDiffCrosstalk = 0;
   Config.RunGeTiming = 0;

   // root output files
   Config.OutPath = "./";
   // Information on energy and wave calibration files
   Config.EnergyCalibrationFile = "./ECal.txt";
   Config.HaveAltEnergyCalibration = 0;
   Config.WaveCalibrationFile = "./WCal.txt";
   Config.HaveWaveCalibration = 0;
   // What to print and how often
   Config.PrintBasic = 1;
   Config.PrintFrequency = PRINT_FREQ;
   Config.PrintVerbose = 0;
   // Limit on number of events
   Config.EventLimit = MAX_EVENTS;
   // ROOT stuff
   Config.ROOT_MaxVirtSize = ROOT_VIRT_SIZE;
   // Where to find the default config file   
   Config.ConfigFile = getenv("GRSISYS");
   Config.ConfigFile += "_Calibrations/Config.txt";

   // Global physics settings
   // ------------------------------------------
   // Properties of waveforms stored in the data
   Config.WaveformSamples = 200;
   Config.WaveInitialSamples = 65;
   Config.WaveFinalSamples = 65;
   // Charge evaluation
   Config.Integration = INTEGRATION;
   // Source information
   float Sources[4][10] = {
      {1173.237, 1332.501},     // 60Co
      {121.7817, 1408.006, 244.6975, 344.2785, 411.116, 778.9040, 964.079, 1112.074, 867.378, 1299.14 },   // 152Eu
      {344.2785, 1408.006, 244.6975, 411.116, 778.9040, 964.079, 1112.074, 867.378, 1299.14},     // 152Eu (no 121)
      {276.398, 356.017, 80.9971, 302.853, 383.851}     // 133Ba
   };
   // Push source lines to vector of vectors.
   // This is horrible.  Need to find another way to initialise these vectors.
   // These values are overwritten if alternatives are given in Config.txt
   vector < float >SourceTemp;
   for (int i = 0; i < 2; i++) {
      SourceTemp.push_back(Sources[SOURCE_60CO][i]);
   }
   Config.Sources.push_back(SourceTemp);
   SourceTemp.clear();
   for (int i = 0; i < 10; i++) {
      SourceTemp.push_back(Sources[SOURCE_152EU_FRONT][i]);
   }
   Config.Sources.push_back(SourceTemp);
   SourceTemp.clear();
   for (int i = 0; i < 9; i++) {
      SourceTemp.push_back(Sources[SOURCE_152EU_BACK][i]);
   }
   Config.Sources.push_back(SourceTemp);
   SourceTemp.clear();
   for (int i = 0; i < 5; i++) {
      SourceTemp.push_back(Sources[SOURCE_133BA][i]);
   }
   Config.Sources.push_back(SourceTemp);
   // Store source names
   Config.SourceNames.push_back("60Co");
   Config.SourceNames.push_back("152Eu");
   Config.SourceNames.push_back("152Eu");
   Config.SourceNames.push_back("133Ba");
   
   // Energy thresholds
   Config.EnergyThresh = 5;     // keV
   Config.ChargeThresh = 100;
   // Array mode
   Config.HighEffMode = 0;

   // Options for Calib() and CalibOffline()
   //---------------------------------------
   // What to do
   Config.CalEnergy = 1;        // Calibrate charge spectra?
   Config.CalWave = 1;          // Calibrate wave spectra?
   Config.Cal2D = 1;            // Create 2D seg-core charge matrices
   Config.CalReport = 1;        // Generate report on fits etc
   Config.CalFile = 0;          // Generate .cal file as used by GRSISpoon 
   memset(&Config.CalList, 1, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
   Config.CalListProvided = 0;
   Config.CalOutputBad = 0;           // output to gain file if calibration is bad?
   Config.CalOverwriteBad = 1;           // Overwrite bad calibration values with mean (useful to get bad seg cals roughly right)
   // Output
   Config.CalOut = "CalibOut.root";     // Name for Tree calibration output file
   Config.CalName = Config.CalOut.substr(0, Config.CalOut.size() - 5);  // name expected for hist calib input. 
      // .root stripped because run number may exist 
   Config.AnaName = "his";
   Config.CalSpecOut = "CalibSpecOut";
   Config.WriteFits = 1;
   // Charge spectra
   Config.ChargeBins = 16384;
   Config.ChargeBins2D = 1024;
   Config.ChargeMax = 1500000;
   Config.WaveChargeMax = 16384;
   // Gain drift/time dependent stuff
   Config.MaxTime = 36000.0;
   Config.TimeBins = 20;
   Config.TimeBinSize = Config.MaxTime / Config.TimeBins;
   Config.FitTempSpectra = 1;
   // Plots
   Config.PlotFits = 0;
   Config.PlotCalib = 0;
   Config.PlotResidual = 0;
   Config.PlotCalibSummary = 1;
   memset(&Config.CalibPlots, 0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
   // Peak Search stuff
   Config.EnSearchThresh = 0.055;
   Config.EnSearchSigma = 10;
   Config.WaveSearchThresh = 0.01;
   Config.WaveSearchSigma = 20;
   // Energy/ch fitting 
   // Estimate of gain for input of fit
   Config.TIGGainEst = 0.1603;  // from average of several TIGRESS HPGe calibrations
   Config.TIGWaveGainEst = 0.6;
   // Extra calibration point at 0 ch = 0 keV
   Config.FitZero = 0;
   Config.ZeroErr = 0.01;
   Config.ForceLinear = 0;         // force linear calibration, even if numlines > 2
   // Peak fitting parameters
   Config.MinFitCounts = 500;       // Minimum counts in whole spectrum for fit to be attempted
   Config.FitWidth_keV = 15;       // Width of region either side of peak to be fitted
   Config.FitBackground = 1;      // 1 = yes, 0 = no.  Should be best to use this all the time but left option there just in case.
   Config.BackWidth_keV = 5;      // Width at each side of fit region to be used as background estimate
   // Initial values for custom fit functions
   Config.EnergySigmaZero = 0.45;
   Config.EnergySigma1MeV = 0.45; 
   Config.WaveSigmaZero   = 1.5; 
   Config.WaveSigma1MeV   = 0.4;
   // Checks on fit quality
   Config.GausHeightMin  = 10.0;
   Config.GausCSPDMax    = 50.0;
   Config.GausSigmaMin   = 250.0;
   // Manual peak select options
   memset(&Config.ManualPeakSelect, 0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
   Config.ManualPeakCorrection = 1;

   // Options for CoincEff()
   //-------------------------
   Config.EffOut = "CoincEffOut.root";
   Config.EffTxtOut = "CoincEffOut.txt";
   Config.EffSimRefFileName = "SimEffRef.txt";
   Config.EffExpRefFileName = "ExpEffRef.txt";
   Config.OutputEff = 1;
   // Plots
   Config.PlotEff = 1;

   // Options for PropXtakl()
   //---------------------------
   Config.PropOut = "PropXTalkOut.root";
   
   // Options for GeTiming
   //------------------------
   // Output
   Config.GeTimingOut = "GeTimingOut.root";
   Config.GeTimingGateCentre = 121.8;
   Config.GeTimingGateWidth = 5.0;
   
   
   return 0;
}

// Override defaults with options from commandd line. 
// Also some options here are required i.e. file and function (--cal,--eff)
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
   int Clover, Crystal, Seg;

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
            cout << endl;
         }
      }
      // Load Configuration file
      // -------------------------------------------
      if (strncmp(argv[i], "-c", 2) == 0) {     // if option is Ecal file
         if (i >= argc - 1 || strncmp(argv[i + 1], "-", 1) == 0) {      // return error if no file
            cout << "No file specified after \"-c\" option." << endl;
            return -1;
         }
         cout << "Checking for configuration file in custom location: " << argv[++i] << "..." << endl;
         if(ReadConfigFile(argv[i]) == 0) {
            ConfigFileLoaded = 1;         // if this is not set default location will also be checked          
            Config.ConfigFile = argv[i];
            cout << "Config Loaded!" << endl << endl;
         }
         else {
            cout << "Failed!" << endl <<endl;  // default will be loaded
         }
         cout << endl;
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
               Config.SourceNumCore.push_back(SOURCE_60CO);
               Config.SourceNumFront.push_back(SOURCE_60CO);
               Config.SourceNumBack.push_back(SOURCE_60CO);
               test = 1;
            }
            if (strncmp(argv[i + 1], "152Eu", 5) == 0 || strncmp(argv[i + 1], "Eu152", 5) == 0 || strncmp(argv[i + 1], "152eu", 5) == 0 || strncmp(argv[i + 1], "eu152", 5) == 0) {     // or is it 152Eu
               Config.SourceNumCore.push_back(SOURCE_152EU_FRONT);
               Config.SourceNumFront.push_back(SOURCE_152EU_FRONT);
               Config.SourceNumBack.push_back(SOURCE_152EU_BACK);
               test = 1;
            }
            if (strncmp(argv[i + 1], "133Ba", 5) == 0 || strncmp(argv[i + 1], "Ba133", 5) == 0 || strncmp(argv[i + 1], "133ba", 5) == 0 || strncmp(argv[i + 1], "ba133", 5) == 0) {     // or is it 133Ba
               Config.SourceNumCore.push_back(SOURCE_133BA);
               Config.SourceNumFront.push_back(SOURCE_133BA);
               Config.SourceNumBack.push_back(SOURCE_133BA);
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
            cout << "Run limited to " << Config.EventLimit << " events." << endl << endl;
         }
      }
      // Alternate output path
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
            // If this option selected, clear detector list
            // Only do this once in case multiple detectors specified
            if (Config.CalListProvided == 0) {
               Config.CalListProvided = 1;
               memset(&Config.CalList, 0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(bool));
            }
            // If keyword cores used..
            if (strncmp(argv[i + 1],"cores",4) == 0 || strncmp(argv[i + 1],"Cores",4) ==0 ) {
               for(Clover=1;Clover<=CLOVERS;Clover++) {
                  for(Crystal=0;Crystal<CRYSTALS;Crystal++) {
                     Config.CalList[Clover-1][Crystal][0] = 1;
                     Config.CalList[Clover-1][Crystal][9] = 1;
                  }
               }
               continue;
            }
            // If keyword segs used...
            if (strncmp(argv[i + 1],"segs",3) == 0 || strncmp(argv[i + 1],"Segs",3) ==0 ) {
               for(Clover=1;Clover<=CLOVERS;Clover++) {
                  for(Crystal=0;Crystal<CRYSTALS;Crystal++) {
                     for(Seg=1;Seg<SEGS+1;Seg++) {
                        Config.CalList[Clover-1][Crystal][Seg] = 1;
                     }
                  }
               }
               continue;
            }
            // else, expect a list of three integers to specify a clover,crystal,seg
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
         // Seg Calibration from core correlation
         if (strncmp(argv[i], "--calsegcore",8) == 0) {
            Config.RunSegCoreCorrelation = 1;
         }
         // calibration
         if (strncmp(argv[i], "--cal", 5) == 0) {
            Config.RunCalibration = 1;
         }
         // efficiency
         if (strncmp(argv[i], "--eff", 5) == 0) {
            Config.RunEfficiency = 1;
         }
         if (strncmp(argv[i], "--speceff", 5) == 0) {
            Config.RunSpecEff = 1;
         }
         // proportional crosstalk
         if (strncmp(argv[i], "--prop", 6) == 0) {
            Config.RunPropCrosstalk = 1;
         }
         if(strncmp(argv[i], "--getim", 7) == 0) {
            Config.RunGeTiming = 1;
         }
      }
   }
   // Look for default configuration file if non given
   if(ConfigFileLoaded == 0) {
      cout << "Checking for configuration file in default location: " << Config.ConfigFile << "..." << endl;
      if(ReadConfigFile(Config.ConfigFile) == 0) {
         cout << "Config Loaded!" << endl << endl;
      }
      else {
         cout << "Failed!" << endl << endl;
         return -1;
      }
      
   }   

   return 0;
}

// Read in data from configuration file.
// Some things e.g. reference values for clover efficienccies will not be used if config file isn't found.
// Other things, like fit parameters, should be set to sensible defaults and overwritten if found in the config file.
int ReadConfigFile(std::string filename) 
{                 
   // Variables
   std::ifstream File;
   std::string Line;
   std::string SubLine;
   int StringPos;
   unsigned int Items = 0;       // count of lines in config file 
   unsigned int Comments = 0;    //   "      comments       "
   unsigned int Other = 0;
   float ValF;
   int ValI;
   int NumRead;
   int Item;
   
   // See if file can be opened.
   File.open(filename,std::ifstream::in);
   if (!File) {  // If not loaded then return with error.
      return 1;
   }
   
   // Looplines and look for known keys .
   while(getline(File,Line)) {
      //cout << Line << endl;
      if(strncmp(Line.c_str(), "#", 1) == 0) {
         Comments += 1;
         continue;
      }
      
      // Reference Values for comparison with measured efficiencies
      // -----------------------------------------------------------
      if (strcmp(Line.c_str(), "SIM_CLOVER_AB_EFF") == 0) {
         if(ReadCloverRef(&File, &Config.Sim_Clover_AB_Eff) > 0) {
            cout << "SIM_CLOVER_AB_EFF reference loaded";
            cout << " (size = " <<  Config.Sim_Clover_AB_Eff.size() << ")"<< endl;
            Items += 1;
            continue;            
         }
      }
      if (strcmp(Line.c_str(), "EXP_CLOVER_AB_EFF") == 0) {
         if(ReadCloverRef(&File, &Config.Exp_Clover_AB_Eff) > 0) {
            cout << "EXP_CLOVER_AB_EFF reference loaded";
            cout << " (size = " <<  Config.Exp_Clover_AB_Eff.size() << ")" << endl;
            Items += 1;
            continue;            
         }
      }
      if (strcmp(Line.c_str(), "SIM_CRYSTAL_EFF") == 0) {
         if(ReadCrystalRef(&File, &Config.Sim_Crystal_Eff) > 0) {
            cout << "SIM_CRYSTAL_EFF reference loaded";
            cout << " (size = " <<  Config.Sim_Crystal_Eff.size() << ")" << endl;
            Items += 1;
            continue;            
         }
      }
      if (strcmp(Line.c_str(), "EXP_CRYSTAL_EFF") == 0) {
         if(ReadCrystalRef(&File, &Config.Exp_Crystal_Eff) > 0) {
            cout << "EXP_CRYSTAL_EFF reference loaded";
            cout << " (size = " <<  Config.Exp_Crystal_Eff.size() << ")" << endl;
            Items += 1;
            continue;            
         }
      }
      
      // Reference values for energy resolutions
      // ---------------------------------------
      if (strcmp(Line.c_str(), "CRYSTAL_FHWM_1332") == 0) {
         if(ReadCrystalRef(&File, &Config.Crystal_FWHM_1332) > 0) {
            cout << "CRYSTAL_FHWM_1332 reference loaded";
            cout << " (size = " <<  Config.Crystal_FWHM_1332.size() << ")" << endl;
            Items += 1;
            continue;            
         }
      }
      
      // Physics/DAQ Settings
      //----------------------
      // Integration used in FPGA charge evaluation
      if (strcmp(Line.c_str(), "INTEGRATION")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%d", &ValI) == 1) {
            Config.Integration = ValI;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      // Source Specification
      //------------------------
      // 60Co
      if (strcmp(Line.c_str(), "60CO")==0) {
         getline(File,Line);
         Config.Sources.at(SOURCE_60CO).clear();  // New data so clear default
         while(sscanf(Line.c_str(),"%f",&ValF)==1) {  // Fails when next item is not a float
            Config.Sources.at(SOURCE_60CO).push_back(ValF);
            StringPos = Line.find(" ");
            Line = Line.substr(StringPos+1,Line.size());
            if(StringPos == -1) {  // No more spaces (this was last item)
               getline(File,Line); // Line pointing to last item, this advances it to first item on next line
               break;
            }
         }
         cout << "60CO Source data loaded (" << Config.Sources.at(SOURCE_60CO).size() << " lines):" << endl;
         for(Item=0;Item<Config.Sources.at(SOURCE_60CO).size();Item++) {
            cout << Config.Sources.at(SOURCE_60CO).at(Item) << " ";
         }
         cout << endl;
      }
      // 152 Eu front segments
      if (strcmp(Line.c_str(), "152EU FRONT")==0) {
         getline(File,Line);
         Config.Sources.at(SOURCE_152EU_FRONT).clear();  // New data so clear default
         while(sscanf(Line.c_str(),"%f",&ValF)==1) {  // Fails when next item is not a float
            Config.Sources.at(SOURCE_152EU_FRONT).push_back(ValF);
            StringPos = Line.find(" ");
            Line = Line.substr(StringPos+1,Line.size());
            if(StringPos == -1) {  // No more spaces (this was last item)
               getline(File,Line); // Line pointing to last item, this advances it to first item on next line
               break;
            }
         }
         cout << "152EU (front) Source data loaded (" << Config.Sources.at(SOURCE_152EU_FRONT).size() << " lines):" << endl;
         for(Item=0;Item<Config.Sources.at(SOURCE_152EU_FRONT).size();Item++) {
            cout << Config.Sources.at(SOURCE_152EU_FRONT).at(Item) << " ";
         }
         cout << endl;
      }
      // 152Eu back segments
      if (strcmp(Line.c_str(), "152EU BACK")==0) {
         getline(File,Line);
         Config.Sources.at(SOURCE_152EU_BACK).clear();  // New data so clear default
         while(sscanf(Line.c_str(),"%f",&ValF)==1) {  // Fails when next item is not a float
            Config.Sources.at(SOURCE_152EU_BACK).push_back(ValF);
            StringPos = Line.find(" ");
            Line = Line.substr(StringPos+1,Line.size());
            if(StringPos == -1) {  // No more spaces (this was last item)
               getline(File,Line); // Line pointing to last item, this advances it to first item on next line
               break;
            }
         }
         cout << "152EU (back) Source data loaded (" << Config.Sources.at(SOURCE_152EU_BACK).size() << " lines):" << endl;
         for(Item=0;Item<Config.Sources.at(SOURCE_152EU_BACK).size();Item++) {
            cout << Config.Sources.at(SOURCE_152EU_BACK).at(Item) << " ";
         }
         cout << endl;
      }
      // 133Ba
      if (strcmp(Line.c_str(), "133BA")==0) {
         getline(File,Line);
         Config.Sources.at(SOURCE_133BA).clear();  // New data so clear default
         while(sscanf(Line.c_str(),"%f",&ValF)==1) {  // Fails when next item is not a float
            Config.Sources.at(SOURCE_133BA).push_back(ValF);
            StringPos = Line.find(" ");
            Line = Line.substr(StringPos+1,Line.size());
            if(StringPos == -1) {  // No more spaces (this was last item)
               getline(File,Line); // Line pointing to last item, this advances it to first item on next line
               break;
            }
         }
         cout << "133Ba Source data loaded (" << Config.Sources.at(SOURCE_133BA).size() << " lines):" << endl;
         for(Item=0;Item<Config.Sources.at(SOURCE_133BA).size();Item++) {
            cout << Config.Sources.at(SOURCE_133BA).at(Item) << " ";
         }
         cout << endl;
      }
      // Calibration/fit parameters
      // ----------------------------
      // Fit temporary spectra during --cal run for gain drift check?
      if (strcmp(Line.c_str(), "FIT_TEMP_SPECTRA")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%d", &ValI) == 1) {
            Config.FitTempSpectra = ValI;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      // Time spectra properties
      if (strcmp(Line.c_str(), "MAX_TIME")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.MaxTime = ValF;
            Config.TimeBinSize = Config.MaxTime / Config.TimeBins;  // Recalculate bin size
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "TIME_BINS")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%d", &ValI) == 1) {
            Config.TimeBins = ValI;
            Config.TimeBinSize = Config.MaxTime / Config.TimeBins;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      // Charge spectra properties
      if (strcmp(Line.c_str(), "CHARGE_BINS")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%d", &ValI) == 1) {
            Config.ChargeBins = ValI;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "CHARGE_BINS2D")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%d", &ValI) == 1) {
            Config.ChargeBins2D = ValI;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "CHARGE_MAX")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.ChargeMax = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "WAVE_CHARGE_MAX")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.WaveChargeMax = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      // Peak Search stuff
      if (strcmp(Line.c_str(), "EN_SEARCH_THRESH")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.EnSearchThresh = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "EN_SEARCH_SIGMA")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.EnSearchSigma = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "WAVE_SEARCH_THRESH")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.WaveSearchThresh = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "WAVE_SEARCH_SIGMA")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.WaveSearchSigma = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "TIG_GAIN_EST")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.TIGGainEst = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "TIG_WAVE_GAIN_EST")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.TIGWaveGainEst = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      // Peak fitting
      if (strcmp(Line.c_str(), "MIN_FIT_COUNTS")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%d", &ValI) == 1) {
            Config.MinFitCounts = ValI;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "FIT_WIDTH_KEV")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.FitWidth_keV = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "BACK_WIDTH_KEV")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.BackWidth_keV = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "FIT_BACKGROUND")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%d", &ValI) == 1) {
            Config.FitBackground = ValI;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "ENERGY_SIGMA_ZERO")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.EnergySigmaZero = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "ENERGY_SIGMA_1MEV")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.EnergySigma1MeV = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "WAVE_SIGMA_ZERO")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.WaveSigmaZero = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "WAVE_SIGMA_1MEV")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.WaveSigma1MeV = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "GAUS_HEIGHT_MIN")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.GausHeightMin = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "GAUS_CSPD_MAX")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.GausCSPDMax = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      if (strcmp(Line.c_str(), "GAUS_SIGMA_MIN")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.GausSigmaMin = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }
      
      // GeTiming Stuff
      //----------------
      if (strcmp(Line.c_str(), "GE_TIMING_GATE_CENTRE")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.GeTimingGateCentre = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      }      
      if (strcmp(Line.c_str(), "GE_TIMING_GATE_WIDTH")==0) {
         getline(File,Line);
         if(sscanf(Line.c_str(), "%f", &ValF) == 1) {
            Config.GeTimingGateWidth = ValF;
            Items += 1;
         }
         else {Other += 1;}
         continue;
      } 
      
      // Other
      //-------
      if (strcmp(Line.c_str(), "")!=0) {
         cout << "WARNING: Unable to idenify configuration item: " << Line.c_str() << " skipping..." << endl;
      }
      Other += 1;
      
   }
   
   cout << "Read " << Items << " items, " << Comments << " comments, " << Other << " others." << endl;
   
   return 0;
}

// Called from ReadConfigFile() to read values by Clover.
// expects format "int float" for clover, Value 
int ReadCloverRef(std::ifstream *File, ReferenceValueMap *Map) {
   
   std::string Line;
   std::vector<int> Vect;
   int Count = 0;
   int Clover;
   float Val;
   while(getline(*File,Line)) {
      if(sscanf(Line.c_str(), "%d %f", &Clover, &Val) == 2) {
         //cout << Line << " done!" << endl;
         Vect.clear();
         Vect.push_back(Clover);
         Map->insert(ReferenceValuePair(Vect,Val));
         Count += 1;
      }
      else {
         return Count;
      }
   }
   return Count;
}

// Called from ReadConfigFile() to read values by Crystal.
// expects format "int int float" for clover, crystal, Value 
int ReadCrystalRef(std::ifstream *File, ReferenceValueMap *Map) {

   std::string Line;
   std::vector<int> Vect;
   int Count = 0;
   int Clover;
   int Crystal;
   float Val;
   while(getline(*File,Line)) {
      if(sscanf(Line.c_str(), "%d %d %f", &Clover, &Crystal, &Val) == 3) {
         //cout << Line << " done!" << endl;
         Vect.clear();
         Vect.push_back(Clover);
         Vect.push_back(Crystal);
         Map->insert(ReferenceValuePair(Vect,Val));
         Count += 1;
      }
      else {
         return Count;
      }
   }
   return Count;
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
