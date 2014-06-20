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

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "Calib.h"
#include "HistCalib.h"
#include "2DHistCalib.h"
#include "Utils.h"

extern TApplication *App;

// globals
TCanvas *cCalib = NULL;


int CorrelateSegsCores() {
   
   // Variables, Constants, etc
   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   std::string Filename;
   TFile *File = NULL;
   char CharBuf[CHAR_BUFFER_SIZE];
   TH1F *Histo = NULL;
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
   
   // Set up TCanvas
   cCalib = new TCanvas("cCalib", "2D Calib", 800, 600);
      
   // Now file should be open, loop Cl,Cr,Seg and load histos
   for(Clover=1; Clover <= CLOVERS; Clover++) {
      for(Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for(Seg=1; Seg<=SEGS; Seg++) { // should looop segs but not cores
            Histo = NULL;
            snprintf(CharBuf,CHAR_BUFFER_SIZE,"TIG%02d%cP%02dx Chg Mat",Clover,Num2Col(Crystal),Seg);
            cout << CharBuf << endl;
            Histo = (TH1F*) File->FindObjectAny(CharBuf);
            if(Histo) {
               cCalib->cd();
               Histo->Draw("colz");
               App->Run(1);
            }
            
         }
      }
   }
   
   
   File->Close();

   return 0;

}











