// g++ CalibOffline.C FitGammaSpectrum.C --std=c++0x -o CalibOffline -O2 `root-config --cflags --libs` -lSpectrum -lgsl -lgslcblas -g

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
#include <TStyle.h>


// TriScope libraries
//#include "TTigFragment.h"

// My libraries
#include "Calib.h"
#include "Main.h"
#include "FitGammaSpectrum.h"


// Sources
float Sources[3][10] = {
   {1173.237, 1332.501},
   {121.7817, 1408.006, 244.6975, 344.2785, 411.116, 778.9040, 964.079, 1112.074},
   {344.2785, 1408.006, 121.7817, 244.6975, 411.116, 778.9040, 964.079, 1112.074}
};

// Fitting stuff
std::string FitOptions = ("RQE");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
std::string Opts;
//float InitialSigma = 0.0;
float InitialSigma = 50.0;
float InitialGain = 0.16;

// Functions called from here:
//int FitSpectrum(TH1F* Histo, SpectrumFit *Fit, int Source , int PlotOn);

TApplication* App;  // Pointer to root environment for plotting etc

TCanvas *c1, *c2, *ctemp;

int main(int argc, char** argv){

   // Variables
   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   int Source = 0;
   int FitSuccess = 0;
   int PlotOn = 0;
   int i;
   char HistName[1024];
   char Charge[10];
   char Colours[] = "BGRW";
   TH1F *GainPlot, *OffsetPlot, *QuadPlot;
   TH1F *GainHist, *OffsetHist, *QuadHist;
   int ItemNum = 0;
   ofstream GainOut;
   ofstream ReportOut;
   
   gStyle->SetOptStat("iouRMen");
   
   // Set up plots
   App= new TApplication("Output", 0, NULL);  // creates root environment for interacting with plots etc
   if(PLOT_FITS || PLOT_CALIB) {
      c1 = new TCanvas();  // Canvas for spectrum plots
      c1->Divide(1,2);
   }
   if(PLOT_CALIB_SUMMARY) {
      c2 = new TCanvas("c2", "Calibration Canvas", 800, 600);  // Canvas for gain plots and histograms
      c2->Divide(2,3);
      c2->Update();   
      ctemp = new TCanvas("ctemp", "Temp Canvas", 100, 50);  // Canvas for gain plots and histograms
         // Having this canvas seems to make c2 division work properly.  I don't know why.  I hate root.
   }
   
   // Files
   if(argc<2) {
      cout << "Filename!" << endl;
      return 0;
   }
   TFile *file = TFile::Open(argv[1]);
   if(OUTPUT_GAIN) {
      GainOut.open("GainsOut.txt");
   }   
   if(OUTPUT_REPORT) {
      ReportOut.open("CalibrationReport.txt");
   }
   
   // Histograms
   GainPlot = new TH1F("Gains","Gain of all fitted channels",1001,-0.5,1000.5);
   GainPlot->GetYaxis()->SetTitle("keV/ch");
   GainPlot->GetXaxis()->SetTitle("Channel");
   OffsetPlot = new TH1F("Offsets","Offset of all fitted channels",1001,-0.5,1000.5);
   OffsetPlot->GetYaxis()->SetTitle("keV");
   OffsetPlot->GetXaxis()->SetTitle("Channel");
   QuadPlot = new TH1F("Quads","Quadratic component of all fitted channels",1001,-0.5,1000.5);
   QuadPlot->GetYaxis()->SetTitle("keV/ch^2");
   QuadPlot->GetXaxis()->SetTitle("Channel");
   
   GainHist = new TH1F("Gain Histogram","Histogram of Gain of all fitted channels",512,0,0.3);
   GainHist->GetXaxis()->SetTitle("keV/ch");
   OffsetHist = new TH1F("Offset Histogram","Histogram of Offset of all fitted channels",512,-5,5);
   OffsetHist->GetXaxis()->SetTitle("keV");
   QuadHist = new TH1F("Quadratic Histogram","Histogram of Quadratic component of all fitted channels",512,-0.000001,0.000001);
   QuadHist->GetXaxis()->SetTitle("keV/ch^2");
   
   if(file){
      cout << argv[1] << " opened." << endl;
   
      for(Clover=0; Clover<CLOVERS; Clover++) {
         for(Crystal=0; Crystal<CRYSTALS; Crystal++) {
            for(Seg=0; Seg<=SEGS+1; Seg++){
               
               switch(Seg) {
                  case 0:
                     sprintf(HistName,"TIG%02d%cN%02da Chg",Clover+1,Colours[Crystal],Seg);
                     Source = SOURCE_NUM_CORE;   
                     break;
                  case 9:
                     sprintf(HistName,"TIG%02d%cN%02db Chg",Clover+1,Colours[Crystal],Seg);
                     Source = SOURCE_NUM_CORE;   
                     break;
                  default:
                     sprintf(HistName,"TIG%02d%cP%02d Chg",Clover+1,Colours[Crystal],Seg);
                     if(Seg<5) {
                        Source = SOURCE_NUM_FRONT;
                     }
                     else {
                        Source = SOURCE_NUM_BACK;
                     }
                     break;
               }              
               
               TH1F *Histo = (TH1F*) file->FindObjectAny(HistName);
               if(Histo) {
                  cout << endl << "------------------------------------" << endl;
                  cout << "Hist " << HistName << " loaded" << endl;
                  
                  SpectrumFit Fit = {};
                  PlotOn = 0;
                  if(PLOT_FITS) {
                     if(PLOT_CLOVER==0 || (Clover+1)==PLOT_CLOVER) {
                        if(PLOT_CRYSTAL==0 || (Crystal+1)==PLOT_CRYSTAL) {
                           if(PLOT_SEG==0 || Seg==PLOT_SEG) {
                              PlotOn=1;
                           }
                        }   
                     }   
                  }
                  
                  // Perform Fit
                  FitSuccess = FitGammaSpectrum(Histo, &Fit, Source, PlotOn);
                  
                  // Increment fit param plots
                  if(FitSuccess > 0) {
                     switch(Crystal) {
                        case 0:
                           ItemNum = (Clover*60)+Seg; 
                           break;
                        case 1:
                           ItemNum = (Clover*60)+20+Seg;
                           break;
                        case 2:
                           ItemNum = (Clover*60)+30+Seg;
                           break;
                        case 3:
                           ItemNum = (Clover*60)+50+Seg;
                           break;
                        default:
                           ItemNum = 1000;
                           break;
                     }                 
                     //cout << "C,C,S = " << Clover << ", " << Crystal << ", " << Seg << "  ItemNum = " << ItemNum << endl;
                     GainPlot->SetBinContent( ItemNum+1, Fit.QuadGainFit[1]);  // ItemNum+1 to skip bin 0 which is underflow
                     OffsetPlot->SetBinContent( ItemNum+1, Fit.QuadGainFit[0]);
                     QuadPlot->SetBinContent( ItemNum+1, Fit.QuadGainFit[2]);   
                     
                     GainHist->Fill(Fit.QuadGainFit[1]);
                     OffsetHist->Fill(Fit.QuadGainFit[0]);
                     QuadHist->Fill(Fit.QuadGainFit[2]);                     
                  }
                  // Now print reports on results of fits and calibration.
                  if(OUTPUT_GAIN) {
                     if(FitSuccess > 0) {
                        GainOut << HistName << " " << Fit.LinGainFit[0] << " +/- " << Fit.dLinGainFit[0];
                        GainOut << " " << Fit.LinGainFit[1] << " +/- " << Fit.dLinGainFit[1] << " " <<   Fit.LinGainFit[2] << endl;
                     }
                     else {
                        GainOut << HistName << " Fail!!!" << endl;
                     }
                  }
                  if(OUTPUT_REPORT) {
                     if(FitSuccess > 0) {
                        // Heading for channel
                        ReportOut << endl << "------------------------------------------" << endl << HistName << endl << "------------------------------------------" << endl << endl;
                        for(i=0; i<NUM_LINES; i++) {
                           if(Fit.FitSuccess[i]==1) {
                              ReportOut << Fit.PeakFits[i].Energy << " keV" << endl;
                              ReportOut << Fit.PeakFits[i].Const << " +/- " <<  Fit.PeakFits[i].dConst << "\t";
                              ReportOut << Fit.PeakFits[i].Mean << " +/- " <<  Fit.PeakFits[i].dMean << "\t";
                              ReportOut << Fit.PeakFits[i].Const << " +/- " <<  Fit.PeakFits[i].dConst << endl;
                              ReportOut << Fit.PeakFits[i].ChiSq << "\t" << Fit.PeakFits[i].NDF << "\t" << Fit.PeakFits[i].ChiSq/Fit.PeakFits[i].NDF;
                              ReportOut << endl << endl;
                           }   
                        }
                        ReportOut << "Linear Solution: Offset = " << Fit.LinGain[0] << " +/- " << Fit.dLinGain[0] << "\t";
                        ReportOut <<  "Gain = " << Fit.LinGain[1] << " +/- " << Fit.dLinGain[1] << endl;
                        ReportOut << "Linear Fit: Offset = " << Fit.LinGainFit[0] << " +/- " << Fit.dLinGainFit[0] << "\t";
                        ReportOut << "Gain = " << Fit.LinGainFit[1] << " +/- " << Fit.dLinGainFit[1] << "\t";
                        ReportOut << "CSPD = " << Fit.LinGainFit[2] << endl;
                        ReportOut << "Quadratic Fit: Offset = " << Fit.QuadGainFit[0] << " +/- " << Fit.dQuadGainFit[0] << "\t";
                        ReportOut << "Gain = " << Fit.QuadGainFit[1] << " +/- " << Fit.dQuadGainFit[1] << "\t";
                        ReportOut << "Quad = " << Fit.QuadGainFit[2] << " +/- " << Fit.dQuadGainFit[2] << "\t";
                        ReportOut << "CSPD = " << Fit.LinGainFit[3] << endl;
                                 
                     }
                     else {
                        ReportOut << endl << "------------------------------------------" << endl << HistName << endl << "------------------------------------------" << endl << endl;
                        ReportOut << "Fail Fail Fail! The calibration has failed!" << endl;
                     }
                  }
               }   
               else {
                  cout << endl << "Hist " << HistName << " failed to load." << endl;
               }
               
            }
         }   
      }      
               
   }
   else {
      cout << "Failed to open " << argv[1] << endl;
   }
   
   if(PLOT_CALIB_SUMMARY) {
      //c2->cd();
      c2->cd(1);
      OffsetPlot->Draw();
      c2->Modified();
      c2->Update();
      //App->Run();
      c2->cd(2);
      OffsetHist->Draw();
      c2->Modified();
      c2->Update();
      c2->cd(3);
      GainPlot->Draw();
      c2->Modified();
      c2->Update();
      c2->cd(4);
      GainHist->Draw();
      c2->Modified();
      c2->Update();
      c2->cd(5);
      QuadPlot->Draw();
      c2->Modified();
      c2->Update();
      c2->cd(6);
      QuadHist->Draw();               
      c2->Modified();
      c2->Update();
      App->Run();
   }
     
   
   if(OUTPUT_GAIN) {
      GainOut.close();
   }   
   if(OUTPUT_REPORT) {
      ReportOut.close();
   }
   
   return 0;
   
}



