// To  compile: g++ SortHistos.C HistCalib.C SegCoreCalib.C Options.C Utils.C -I$GRSISYS/include --std=c++0x -o SortHistos  -O0 `root-config --cflags --libs`  -lSpectrum -g
using namespace std;
// C/C++ libraries:
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <time.h>

// ROOT libraries:
#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TCutG.h>
#include <TTreeIndex.h>
#include <TTreePlayer.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TStyle.h>

// TriScope libraries
//#include "TTigFragment.h"

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "Calib.h"
#include "HistCalib.h"
#include "SegCoreCalib.h"
#include "Utils.h"

TApplication *App;              // Pointer to root environment for plotting etc

int main(int argc, char **argv)
{
   int i;

   // Timing
   TStopwatch StopWatch;
   StopWatch.Start();

   // Set default and read custom options
   LoadDefaultSettings();
   if (ReadCommandLineSettings(argc, argv) < 0) {
      cout << "Failed to configure the run - exiting!" << endl;
      return -1;
   }
   
   // Load any alternate calibration information 
   if (Config.HaveAltEnergyCalibration) {
      int NumCal;
      NumCal = ReadCalibrationFile(Config.EnergyCalibrationFile, &Config.EnCalibNames, &Config.EnCalibValues);
      if (Config.PrintBasic) {
         cout << "Alternate energy calibratrion values: " << endl;
         for (i = 0; i < NumCal; i++) {
            cout << i;
            cout << ": " << Config.EnCalibNames.at(i);
            cout << " " << Config.EnCalibValues.at(i)[0] << " " << Config.EnCalibValues[i][1] << " " << Config.EnCalibValues[i][2] << endl;
         }
      }
   }
   // Load gain coefficients for waveforms
   if (Config.HaveWaveCalibration) {
      int NumCal;
      NumCal = ReadCalibrationFile(Config.WaveCalibrationFile, &Config.WaveCalibNames, &Config.WaveCalibValues);
      if (Config.PrintBasic) {
         cout << "Wave energy calibratrion values: " << endl;
         for (i = 0; i < NumCal; i++) {
            cout << i;
            cout << ": " << Config.WaveCalibNames.at(i);
            cout << " " << Config.WaveCalibValues.at(i)[0] << " " << Config.WaveCalibValues[i][1] << " " << Config.WaveCalibValues[i][2] <<
                endl;
         }
      }
   }
   
   // Set options for histo stats
   gStyle->SetOptStat("iouRMen");

   // create root environment for interacting with plots etc
   App = new TApplication("Output", 0, NULL);

   // Check what we are supposed to be doing, call function
   if (Config.RunSpecCal == 1) {
      if (CalibrateFiles() > 0) {
         cout << "CalibrateFiles() failed." << endl;
         return 1;
      }
   }

   if(Config.RunSegCoreCorrelation == 1) {
      if (!Config.HaveAltEnergyCalibration) {
         cout << "Calibration file for cores required for seg-core correlation calibration of seg energies." << endl;
      }
      if (!Config.HaveWaveCalibration) {
         cout << "Wave calibration file for cores required for seg-core correlation calibration of seg wave energies." << endl;
      }
      if (!Config.HaveAltEnergyCalibration && !Config.HaveWaveCalibration) {
         cout << "Error - Calibration files required for Seg-Core correlation calibration." << endl;
         return 1;
      }
      if (SegCoreCalib() > 0) {
         cout << "CorrelateCoresSegs() failed." << endl;
         return 1; 
      }
   }

   if (Config.RunSpecEff == 1) {
      cout << "Relative efficiency fitting not yet implemented." << endl;
   }
   
   if(Config.PrintBasic) {
      cout << "SortHistos completed in " << StopWatch.RealTime() << " seconds." << endl;
   }
   
   return 0;
}
