// To  compile: g++ SortHistos.C HistCalib.C Options.C Utils.C -I$GRSISYS/include --std=c++0x -o SortHistos  -O0 `root-config --cflags --libs`  -lSpectrum -g
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
#include "2DHistCalib.h"
#include "Utils.h"

TApplication *App;              // Pointer to root environment for plotting etc

int main(int argc, char **argv)
{
   // Timing
   TStopwatch StopWatch;
   StopWatch.Start();

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

   // Check what we are supposed to be doing, call function
   if (Config.RunSpecCal == 1) {
      if (CalibrateFiles() > 0) {
         cout << "CalibrateFiles() failed." << endl;
         return 1;
      }
   }

   if(Config.RunSegCoreCorrelation == 1) {
      if (CorrelateSegsCores()> 0) {
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
