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
#include <TProfile.h>

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "Calib.h"
#include "HistCalib.h"
#include "SegCoreCalib.h"
#include "Utils.h"

extern TApplication *App;

// globals
TCanvas *cCalib = NULL;


int SegCoreCalib() {
   
   // Variables, Constants, etc
   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   std::string Filename;
   TFile *File = NULL;
   std::string tempstring;
   char CharBuf[CHAR_BUFFER_SIZE];
   TH2F *Histo = NULL;
   ofstream SegCoreCalOut;
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
   
   // Output file
   tempstring = Config.OutPath + "SegCoreCalOut.txt";
   SegCoreCalOut.open(tempstring.c_str());
   
   // Set up TCanvas
   cCalib = new TCanvas("cCalib", "2D Calib", 800, 600);
      
   // Now file should be open, loop Cl,Cr,Seg and load histos
   for(Clover=1; Clover <= CLOVERS; Clover++) {
      for(Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for(Seg=1; Seg<=SEGS; Seg++) { // should looop segs but not cores
            Histo = NULL;
            snprintf(CharBuf,CHAR_BUFFER_SIZE,"TIG%02d%cP%02dx Chg Mat",Clover,Num2Col(Crystal),Seg);
            cout << CharBuf << endl;
            Histo = (TH2F*) File->FindObjectAny(CharBuf);
            if(Histo) {
               cCalib->cd();
               
               // -------------------------------------------------------------------------
               // Now to fit matrices and extract transform from core to seg calibration.
               // Two strategies for this seem apparant from the ROOT documentation:
               //    - Create TProfile of matrix with Histo->ProfileX() then fit that as 1D
               //    - Draw matrix and then extract TGraph from that with something like:
               //       TGraph *Graph = (TGraph*) gPad->GetPrimitive("Graph");
               // So far I have had more look with the first method.  Either way there
               // are lots of counts (~10%) with SegCharge<CoreCharge so need to get rid
               // of those first.  Will try using a threshold.
               // -------------------------------------------------------------------------               
               
               Histo->Draw();
               
               App->Run(1);
               
               // Subtract background
               int x,y;
               float Bgnd = 10.0;
               float Val = 0.0;
               for(x=0;x<Histo->GetNbinsX();x++) {
                  for(y=0;y<Histo->GetNbinsY();y++) {
                     Val = Histo->GetBinContent(x,y);
                     if(Val>Bgnd) {
                        Histo->SetBinContent(x,y,Val-Bgnd);
                     }
                     else {
                        Histo->SetBinContent(x,y,0.0);
                     }
                  }
               }
               
               Histo->Draw();
               
               App->Run(1);
               
               TProfile *profx = Histo->ProfileX();
               
               profx->Fit("pol2");
               
               
               SegCoreCalOut << CharBuf << " " << endl;
               /*for(Param=0;Param<2;Param++){
                  SegCoreCalOut << 
               }*/
               
               
               profx->Draw();
               
               App->Run(1);
               
            }
            
         }
      }
   }
   
   SegCoreCalOut.close();
   File->Close();

   return 0;

}











