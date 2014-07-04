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
   int Param;
   std::string Filename;
   TFile *File = NULL;
   std::string tempstring;
   char CharBuf[CHAR_BUFFER_SIZE];
   std::string CoreName;
   std::string SegName;
   TH2F *Histo = NULL;
   ofstream SegCoreCalOut;
   // Fitting stuff
   std::string FitOptions = ("RQE");
      // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
   // Calibration
   std::vector<float> Coeffs;      // Seg-core correlation Coeffs
   std::vector<float> CoreCoeffs;  // Core Calibration Coeffs
   std::vector<float> SegCoeffs;   // Seg calibration coefficiets
   
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
            snprintf(CharBuf,CHAR_BUFFER_SIZE,"TIG%02d%cN00a",Clover,Num2Col(Crystal));
            CoreName = CharBuf;
            snprintf(CharBuf,CHAR_BUFFER_SIZE,"TIG%02d%cP%02dx Chg Mat",Clover,Num2Col(Crystal),Seg);
            cout << CharBuf << endl;
            SegName = CharBuf;
            
            Histo = (TH2F*) File->FindObjectAny(SegName.c_str());
            if(Histo) {
               
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
               
               if(Config.PlotSegCoreCal == 1) {
                  cCalib->cd();
                  Histo->Draw("colz");
                  App->Run(1);
               }
               
               // Variables needed here
               int x,y;
               float Bgnd;
               float Val;
               float Int;
               int CalChan;
               bool NewCoeffFound;
               int i;
               TF1 *ProfileFit;
               float Min, Max;
               
               // skip if stats too low
               cout << "Counts: " << Histo->Integral() << endl;
               if(Histo->Integral() < Config.MinFitCounts) {
                  //continue;
               }
               
               // Subtract background to remove values with Eseg < Ecore
               Bgnd = 10.0; // If sticking with this method, this value should be found from matrix not hard coded.
               Val = 0.0;
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
               
               
               if(Config.PlotSegCoreCal == 1) {
                  cCalib->cd();
                  Histo->Draw("colz");
                  App->Run(1);
               }
               
               TProfile *ProfX = Histo->ProfileX();
               
               Min = ProfX->GetMinimum();
               Max = ProfX->GetMaximum();
               
               if(Config.PlotSegCoreCal == 1) {
                  cCalib->cd();
                  ProfX->Draw();
                  App->Run(1);
               }
               
               // Set function based on required order of Seg-Core correlation
               switch(Config.SegCoreFitOrder) {
               case 0:
                  ProfileFit = new TF1("Gain", "x*[0]", Min, Max);
                  break;
               case 1:
                  ProfileFit = new TF1("Pol1", "[0]+(x*[1])", Min, Max);
                  break;
               case 2:
                  ProfileFit = new TF1("Pol2", "[0]+(x*[1])+(x*x*[2])", Min, Max);
                  break;
               }
                 
               
               bool FitTest = ProfX->Fit(ProfileFit,FitOptions.c_str());
               
               // Test fit success
               if(FitTest>0) {
                  if(Config.PrintBasic) {
                     cout << "Fit of TProfile failed." << endl;
                  }
                  SegCoreCalOut << SegName.c_str() << " Fit of TProfile failed!" << endl;
                  continue;
               }
               
               // Now generate output of correlation fit        
               Coeffs.clear();
               SegCoreCalOut << SegName.c_str() << " ";
               for(Param=0;Param<=Config.SegCoreFitOrder;Param++){
                  SegCoreCalOut << ProfileFit->GetParameter(Param) << " ";
                  Coeffs.push_back(ProfileFit->GetParameter(Param));
               }

               // Now output effective seg gain coefficients
               // ------------------------------------------
               CoreCoeffs.clear();
               SegCoeffs.clear();
               
               // Loop calibration  core coeffs
               NewCoeffFound = 0;
               for (CalChan = 0; CalChan < Config.EnCalibNames.size(); CalChan++) {
                  //cout << "CN: " << CoreName << " CoeffN: " << Config.EnCalibNames[CalChan].c_str() << endl;
                  if (strncmp(Config.EnCalibNames[CalChan].c_str(), CoreName.c_str(), 9) == 0) {   // bug!  this will match the first core
                     // name it finds to either a OR b.  Compare 10 chars woud work but then case sensitivity isses on the x/a/b 
                     // at the end.  Don't really need second core energy right now so I will come back to this later
                     NewCoeffFound = 1;
                     break;
                  }
               }
               if (NewCoeffFound == 1) {        // If a new set of coeffs was found, then calibrate
                  CoreCoeffs = Config.EnCalibValues.at(CalChan);
               } else {         // else use the existing calibration
                  if(Config.PrintBasic) {
                     cout << "Calibration of " << SegName << " failed as no core coeffs found." << endl;
                     continue;
                  }
               }
               
               // Check there are 3 core coeffs 
               if(CoreCoeffs.size() != 3) {
                  if(CoreCoeffs.size() == 2) {
                     CoreCoeffs.push_back(0.0);
                  }
                  else {
                     if(Config.PrintBasic) {
                        cout << "Calibration of " << SegName << " failed. Wrong number of core coeffs." << endl;
                        continue;
                     }
                  }               
               }
               
               Int = Config.Integration / Config.Dispersion;
               
               switch(Config.SegCoreFitOrder) {
               case 0:  // Qcore = k * Qseg
                  // s0 = c0
                  Val = CoreCoeffs.at(0);
                  SegCoeffs.push_back(Val);
                  // s1 = c1 * k
                  Val = CoreCoeffs.at(1) * Coeffs.at(0);
                  SegCoeffs.push_back(Val);
                  // s2 = c2 * k**2
                  Val = CoreCoeffs.at(2) * pow(Coeffs.at(0),2);
                  SegCoeffs.push_back(Val);
                  break;
               case 1:  // Qcore = k0 + k1 * Qseg
                  // s0 = c0 + (C1 * k0 / int) + (c2 * K0**2  / int**2)
                  Val = CoreCoeffs.at(0) + (CoreCoeffs.at(1) * Coeffs.at(0) / Int);
                  Val += (CoreCoeffs.at(2) * pow(Coeffs.at(0),2) / pow(Int,2));        
                  SegCoeffs.push_back(Val);
                  // s1 = (c1 * k1) + (2 * c2 * k0 * k1 / int)   
                  Val = (CoreCoeffs.at(1) * Coeffs.at(1));
                  Val += (2*Coeffs.at(0)*Coeffs.at(1)*CoreCoeffs.at(2) / Int);  
                  SegCoeffs.push_back(Val);
                  // s2 = c2 * k1**2
                  Val = CoreCoeffs.at(2) * pow(Coeffs.at(1),2);
                  SegCoeffs.push_back(Val);
                  break;
               case 2:  // Qcore = k0 + (k1*Qseg) + (k2 * Qseg**2)
                  // s0 = c0 + (c1 * k0 /int) + (c2 * K0**2  / int**2)
                  Val = CoreCoeffs.at(0) + (CoreCoeffs.at(1) * Coeffs.at(0) / Int);
                  Val += (CoreCoeffs.at(2) * pow(Coeffs.at(0),2) / pow(Int,2)); 
                  SegCoeffs.push_back(Val);
                  // s1 = (c1*k1) + (2 * c2 * k0 * k1 / int) 
                  Val = (CoreCoeffs.at(1) * Coeffs.at(1));
                  Val += (2.0*Coeffs.at(0)*Coeffs.at(1)*CoreCoeffs.at(2) / Int);  
                  SegCoeffs.push_back(Val);
                  // s2 = (c2 * k1**2) + (c1 * k2 * int) + (2 * c2 * k0 * k2)
                  Val = CoreCoeffs.at(2) * pow(Coeffs.at(1),2);
                  Val += CoreCoeffs.at(1) * Coeffs.at(2) * Int;
                  Val += 2.0 * CoreCoeffs.at(2) * Coeffs.at(0) * Coeffs.at(2);
                  SegCoeffs.push_back(Val);
                  // s3 = 2 * c2 * k1 * k2 * int
                  Val = 2.0 * CoreCoeffs.at(2) * Coeffs.at(1) * Coeffs.at(2) * Int;
                  SegCoeffs.push_back(Val);
                  // s4 = c2 * k2**2 * int**2
                  Val = CoreCoeffs.at(2) * pow(Coeffs.at(2),2) * pow(Int,2);
                  SegCoeffs.push_back(Val);
                  break;
               }
               
               for(Param=0;Param<SegCoeffs.size();Param++){
                  SegCoreCalOut << SegCoeffs.at(Param) << " ";
               }
               SegCoreCalOut << endl;
               
               if(Config.PlotSegCoreCal == 1) {
                  cCalib->cd();
                  ProfX->Draw();
                  App->Run(1);
               }
               
            }
            
         }
      }
   }
   
   SegCoreCalOut.close();
   File->Close();

   return 0;

}











