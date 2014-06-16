// C/C++ libraries:
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
using namespace std;
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
//#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TFolder.h>
#include <TRandom3.h>

// TriScope libraries
//#include "TTigFragment.h"

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "Calib.h"
#include "HistCalib.h"
#include "2DHistCalib.h"
#include "Utils.h"


int CorrelateSegsCores() {
   
   // Variables, Constants, etc
   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   std::string Filename;
   TFile *File = NULL;
   
   // Check inputs/configuration
   if(Config.files.size() != 1) {
      cout << "One file expected for seg calibration by seg-core correlation." << endl;
      return 1;
   }

   Filename = Config.files.at(0);
   
   // Input File
   File = TFile::Open(Filename.c_str(), "READ");
   if (File->IsOpen()) {
      if (Config.PrintBasic) {
         cout << Filename << " opened!" << endl;
      }
   } else {
      if (Config.PrintBasic) {
         cout << "Failed to open " << Filename << "!" << endl;
      }
      return 1;
   }
   
      
   
   File->Close();

   return 0;

}











