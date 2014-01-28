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


// TriScope libraries
#include "TTigFragment.h"
//#include "TFSPC_Info.h"
//#include "TSharc.h"
//#include "TTigress.h"
//#include "TRf.h"
//#include "TTriFoil.h"

// My libraries
#include "Calib.h"
#include "Main.h"
//#include "FitSpectrum.h"

// File pointers:
static TFile *outfile = 0;
static TCanvas *c1, *c2, *ctemp;

extern TApplication* App;

// Spectra Pointers:
static TH1F *hCharge[CLOVERS][CRYSTALS][SEGS+1] = {};  // primary Cores and segs
static TH1F *hMidasTime = 0;
static TH1F *hCrystalChargeTemp[CRYSTALS * CLOVERS] = {};  // Only doing these guys for the cores
static TH1F *hCrystalGain[CRYSTALS * CLOVERS] = {};    
static TH1F *hCrystalOffset[CRYSTALS * CLOVERS] = {};



// Structure for fits and gain for one spectrum
struct SpectrumFit   {
   FitResult PeakFits[NUM_LINES];
   int FitSuccess[NUM_LINES];
   float LinGain[2];  // [O,G]
   float dLinGain[2];
   float LinGainFit[3]; // [O,G,CSPD]
   float dLinGainFit[2]; // [dO,dG]
   float QuadGainFit[4]; // [O,G,Q,CSPD]
   float dQuadGainFit[4]; // [dO,dG,dQ]
};



// Sources
float Sources[3][10] = {
   {1173.237, 1332.501},
   {121.7817, 1408.006, 244.6975, 344.2785, 778.9040, 964.079},
   {344.2785, 1408.006, 121.7817, 244.6975, 778.9040, 964.079}
};

// Fitting stuff
std::string FitOptions = ("RQE");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
std::string Opts;
//float InitialSigma = 0.0;
float InitialSigma = 50.0;
float InitialGain = 0.16;


// Functions called from main:   
void InitCalib();
void Calib(std::vector<TTigFragment> &ev);
void FinalCalib();
// Functions called from here:
void ResetTempSpectra();
//int FitSpectrum(TH1F* Histo, float* Gain, float* Offset, float* dGain, float* dOffset);
int FitSpectrum(TH1F* Histo, SpectrumFit *Fit, int Source , int PlotOn);

void InitCalib() {
 
   char Colours[] = "BGRW";

    // Initialise output file                                       
   outfile = new TFile("CalibOut.root","RECREATE");

   char name[512], title[512];
   int Clover, Crystal, Seg, NumBins;


   for(Clover=0;Clover<CLOVERS;Clover++) {
      for(Crystal=0;Crystal<CRYSTALS;Crystal++) {
         
         sprintf(name,"TIG%02d%c Tmp Chg",Clover+1,Colours[Crystal]);
         sprintf(title,"TIG%02d%c Temp Core Charge (arb)",Clover+1,Colours[Crystal]);
         hCrystalChargeTemp[((Clover*4)+Crystal)] = new TH1F(name,title,CHARGE_BINS,0,CHARGE_MAX);
         sprintf(name,"TIG%02d%c Gain",Clover+1,Colours[Crystal]);
         sprintf(title,"TIG%02d%c Gain vs Time",Clover+1,Colours[Crystal]);
         hCrystalGain[((Clover*4)+Crystal)]   = new TH1F(name,title,TIME_BINS,0,MAX_TIME);
         sprintf(name,"TIG%02d%c Offset",Clover+1,Colours[Crystal]);
         sprintf(title,"TIG%02d%c Offset vs Time",Clover+1,Colours[Crystal]);
         hCrystalOffset[((Clover*4)+Crystal)]   = new TH1F(name,title,TIME_BINS,0,MAX_TIME);
         
         for(Seg = 0; Seg<=SEGS; Seg++) {
            //cout << "Clov,Col,Seg :" << Clover+1 << ", " << Colours[Crystal] << ", " << Seg << endl;
            if(Seg==0) {
               sprintf(name,"TIG%02d%cN%02d Chg",Clover+1,Colours[Crystal],Seg);
               sprintf(title,"TIG%02d%cN%02d Charge (arb)",Clover+1,Colours[Crystal],Seg);
            }
            
            else {
               
               sprintf(name,"TIG%02d%cP%02d Chg",Clover+1,Colours[Crystal],Seg);
               sprintf(title,"TIG%02d%cP%02d Charge (arb)",Clover+1,Colours[Crystal],Seg);
            }
            hCharge[Clover][Crystal][Seg] = new TH1F(name,title,CHARGE_BINS,0,CHARGE_MAX);
            
         }
      }
   }
   sprintf(name,"Midas Time");
   sprintf(title,"Midas Timestamps (s)");
   hMidasTime = new TH1F(name,title,TIME_BINS,0,MAX_TIME);

   cout << "Searching for core Peaks: " << Sources[SOURCE_NUM_CORE][0] << "kev and " << Sources[SOURCE_NUM_CORE][1] << "keV (Ratio " << Sources[SOURCE_NUM_CORE][0]/Sources[SOURCE_NUM_CORE][1] << ")" << endl;
   cout << "Searching for front seg peaks: " << Sources[SOURCE_NUM_FRONT][0] << "kev and " << Sources[SOURCE_NUM_FRONT][1] << "keV (Ratio " << Sources[SOURCE_NUM_FRONT][0]/Sources[SOURCE_NUM_FRONT][1] << ")" << endl;
   cout << "Searching for back seg peaks: " << Sources[SOURCE_NUM_BACK][0] << "kev and " << Sources[SOURCE_NUM_BACK][1] << "keV (Ratio " << Sources[SOURCE_NUM_BACK][0]/Sources[SOURCE_NUM_BACK][1] << ")" << endl;
         

   if(PLOT_FITS || PLOT_CALIB) {
      c1 = new TCanvas();  // Canvas for spectrum plots
      c1->Divide(1,2);
   }
   
   //if(PLOT_CALIB) {
     // c2 = new TCanvas();
   //}
   
}


void Calib(std::vector<TTigFragment> &ev) {

   //cout << "HERE!!!!!! " << ev.size() << endl;

   //Variables
   int i, j, k;
   int Slave, Port, Chan;
   int Crystal, Clover, SpecNum;
   int FitSuccess = 0;
   std::string name;
   
   time_t MidasTime;  
   static time_t StartTime;
   static int FirstEvent=1;
   double RunTimeElapsed;
   static double FitTimeElapsed = 0.0;
   
   int TimeBin = 0;
   float TB = 0.0;
   
   int PlotOn = 0;
   int Source = 0;
   
   for(i=0;i<ev.size();i++) {
   
      //Get time of first fragment
      if(FirstEvent==1) {
         StartTime = ev[i].MidasTimeStamp; //ev[i].MidasTimeStamp;
         FirstEvent = 0;
         cout << "MIDAS time of first fragment: " << ctime(&StartTime) << endl;
      }
      // Insert test for time earlier than StartTime.....
      if(ev[i].MidasTimeStamp < StartTime) {
         StartTime = ev[i].MidasTimeStamp; //ev[i].MidasTimeStamp;
         cout << "Earlier event found! Updating time of first fragment to: " << ctime(&StartTime) << endl;
      }
   
      Slave = ((ev[i].ChannelAddress & 0x00F00000) >> 20);
      Port  = ((ev[i].ChannelAddress & 0x00000F00) >> 8);
      Chan  =  (ev[i].ChannelAddress & 0x000000FF);     
      name	=	ev[i].ChannelName;
      //cout << "Slave, Port, Chan = " << Slave << ", " << Port << ", " << Chan << "\t" << name << endl;
      
      Mnemonic mnemonic;
      if(name.size() >= 10) {ParseMnemonic(&name,&mnemonic);}
      else {
         cout << "Bad mnemonic size: This shouldn't happen if the odb is correctly configured!" << endl;
         continue;   
      }
   
      // If TIGRESS then fill energy spectra
      if(mnemonic.system == "TI") {
         // Determine Crystal
         char Colour = mnemonic.arraysubposition.c_str()[0];
         Crystal = Col2Num(Colour);
         if(Crystal == -1) {
            cout << "Bad Colour: " << Colour << endl;
            continue;
         }

         // Determine Clover position
         Clover = mnemonic.arrayposition;
         
         // If Core
         if(mnemonic.segment == 0) {
            if(mnemonic.outputsensor.compare("a") || mnemonic.outputsensor.compare("A")) {  // if this is the primary core output
            // Increment raw charge spectra
            hCharge[Clover-1][Crystal][0]->Fill(ev[i].Charge);
            hCrystalChargeTemp[(((Clover-1)*4)+Crystal)]->Fill(ev[i].Charge);
            }
         }
         else{
            if(mnemonic.segment < 9) {
               //cout << "Clov, Crys, Seg, Chg: " << Clover << " ," << Crystal << ", " << mnemonic.segment << ", " << ev[i].Charge << endl;
               if(ev[i].Charge > 0) {
                  hCharge[Clover-1][Crystal][mnemonic.segment]->Fill(ev[i].Charge);  // Fill segment spectra
               }
            }
         }
         
         // Now Get time elapsed in run
         MidasTime = ev[i].MidasTimeStamp;
         //cout << "Time: " << ctime(&MidasTime) << endl;
         RunTimeElapsed = difftime(MidasTime,StartTime);
         hMidasTime->Fill(RunTimeElapsed);
         //cout << "Time: " << ctime(&MidasTime) << endl;
         // Have we moved on to a new time period?
         if((RunTimeElapsed - FitTimeElapsed >= TIME_BIN_SIZE) && FIT_TEMP_SPECTRA) {
         
            TB = ((MidasTime - StartTime)/MAX_TIME);
            TimeBin = TB*TIME_BINS;
            
            cout << "MT: " << MidasTime << " ST: " << StartTime << " MT: " << MAX_TIME << " TIME_BINS: " << TIME_BINS << " TB: " << TB << " TimeBin: " << TimeBin << endl;
            
            if(NOTIFY_TIME_BIN) {
               cout << "Start Of Run:      " << ctime(&StartTime);
               cout << "Start of this fit: " << FitTimeElapsed << endl;
               cout << "Time of this frag: " << ctime(&MidasTime);
               cout << "RunTimeElapsed:" << RunTimeElapsed << endl;
               cout << "TimeBin: " <<  TimeBin << endl;
            }
            FitTimeElapsed = RunTimeElapsed; 
        

            // Loop crystals, fit spectra core , write gains and offsets to spectrum
            for(j=0; j<CLOVERS; j++) {
               for(k=0; k<CRYSTALS; k++) {
                  if(VERBOSE){
                     cout << "Clov: " << j << " Crys: " << k;
                  }
                  
                  SpectrumFit Fit = {};
                  PlotOn = 0;
                  Source = SOURCE_NUM_CORE;
                  //FitSuccess = FitSpectrum(hCrystalChargeTemp[(j*4)+k], &Gains[(j*4)+k], &Offsets[(j*4)+k], &dGains[(j*4)+k], &dOffsets[(j*4)+k]);
                  FitSuccess = FitSpectrum(hCrystalChargeTemp[(j*4)+k], &Fit, Source, PlotOn);
                  
                  if(FitSuccess > -1) {
                     hCrystalGain[(j*4)+k]->SetBinContent(TimeBin,Fit.LinGainFit[1]);
                     hCrystalOffset[(j*4)+k]->SetBinContent(TimeBin,Fit.LinGainFit[0]);
                  }   
                  else{
                     if(VERBOSE) {
                        switch(FitSuccess) {
                           case -1: continue;
                           case -2: continue;
                           case -3: continue;
                           case -4: continue;
                           default: continue;
                           break;                        
                        }
                     }
                  }
               }
            }
            // Reset temp spectra
            ResetTempSpectra();
         }
      }  
   }  
}



void FinalCalib() {

   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   int FitSuccess = 0;
   int PlotOn = 0;
   int Source = 0;
   int i;
   ofstream GainOut;
   ofstream ReportOut;
   char HistName[1024];
   
   TH1F *GainPlot, *OffsetPlot, *QuadPlot;
   TH1F *GainHist, *OffsetHist, *QuadHist;
   int ItemNum;
   
   if(PLOT_CALIB_SUMMARY) {
      c2 = new TCanvas("c2", "Calibration Canvas", 800, 600);  // Canvas for gain plots and histograms
      c2->Divide(2,3);
      c2->Update();   
      ctemp = new TCanvas("ctemp", "Temp Canvas", 100, 50);  // Canvas for gain plots and histograms
         // Having this canvas seems to make c2 division work properly.  I don't know why.  I hate root.
   }

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
   
   // Write spectra to file
   outfile->cd();
   for(Clover=0;Clover<CLOVERS;Clover++) {
      for(Crystal=0;Crystal<CRYSTALS;Crystal++) {
         for(Seg=0; Seg<=SEGS; Seg++) {
            hCharge[Clover][Crystal][Seg]->Write();
         }
         hCrystalChargeTemp[((Clover*4)+Crystal)]->Write();
         hCrystalGain[((Clover*4)+Crystal)]->Write();
         hCrystalOffset[((Clover*4)+Crystal)]->Write();
      }
   }
   hMidasTime->Write();
   
   // Files for gain/fit output
   if(OUTPUT_GAIN) {
      GainOut.open("GainsOut.txt");
   }   
   if(OUTPUT_REPORT) {
      ReportOut.open("CalibrationReport.txt");
   }
   
   // Fit full-run energy spectra
   if(FIT_FINAL_SPECTRA) {
      if(VERBOSE) {
         cout << endl << "Now fitting energy spectra from whole run..." << endl;
      }
       
   
      for(Clover=0; Clover<CLOVERS; Clover++) {
         for(Crystal=0; Crystal<CRYSTALS; Crystal++) {
            for(Seg=0; Seg<=SEGS; Seg++){
            
               if(VERBOSE){
                  cout << endl << "--------------------------------------" << endl;
                  cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                  cout << "--------------------------------------" << endl;
               }
               
               SpectrumFit Fit = {};
               PlotOn = 0;
               if(Seg==0) {
                  Source = SOURCE_NUM_CORE;
                  sprintf(HistName,"TIG%02d%cN%02d Chg",Clover+1,Num2Col(Crystal),Seg);
               }
               else {
                  sprintf(HistName,"TIG%02d%cP%02d Chg",Clover+1,Num2Col(Crystal),Seg);
                  if(Seg < 5) {
                     Source = SOURCE_NUM_FRONT;
                  }
                  else {
                     Source = SOURCE_NUM_BACK;
                  }
               }
               
               // Perform Fit
               FitSuccess = FitSpectrum(hCharge[Clover][Crystal][Seg], &Fit, Source, PlotOn);
               
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
         }
      }
   
      if(OUTPUT_GAIN) {
         GainOut.close(); 
      }
      if(OUTPUT_REPORT) {
         ReportOut.close();
      }
      
      
   }
  
   if(PLOT_CALIB_SUMMARY) {
      c2->cd(1);
      OffsetPlot->Draw();
      c2->Modified();
      c2->Update();
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
   
   // Write gain histograms to file
   outfile->cd();
   GainPlot->Write();
   OffsetPlot->Write();
   QuadPlot->Write();
   GainHist->Write();
   OffsetHist->Write();
   QuadHist->Write();   
   
   outfile->Close();
   
}


// 

void ResetTempSpectra() {
   int Clover,Crystal;
   for(Clover=0;Clover<CLOVERS;Clover++) {
      for(Crystal=0;Crystal<CRYSTALS;Crystal++) {
         hCrystalChargeTemp[((Clover*4)+Crystal)]->Reset();
      }
   }
}

int FitSpectrum(TH1F* Histo, SpectrumFit *Fit, int Source, int PlotOn) {

   int Peak, Peak1, Peak2, NumPeaks, BestPeak1, BestPeak2, PeakFound;  // for looping peaks and finding correct ones
   float Ratio, Diff, BestDiff;  // test quality of peak match
   float Min, Max, BackEst;
   int i;
   int MinBin, MaxBin, BackBins;
   Float_t* PeakPositions;  // for output of root peak search
   int Line, LinesUsed;   
   float FitCentre, G, O, dG, dO;  // Fits and values for fast 2 point calibration.
   float Energies[NUM_LINES], dEnergies[NUM_LINES], Centroids[NUM_LINES], dCentroids[NUM_LINES]; // Data for full calibration
   float En1, En2;   // line energies for two point calib    
   En1 = Sources[Source][0];
   En2 = Sources[Source][1];
   float IdealRatio = En1/En2;  
   float Chg1, Chg2, dChg1, dChg2;  // charge and erros for 2 point         
   FitResult FitRes[NUM_LINES];  // store full fit results
   int Integral;  // counts in spectrum
   TSpectrum *Spec = new TSpectrum();
   
   if(Histo){
      Integral = Histo->Integral();
   }
   else {
      cout << "Bad histo! " << endl;
   }
   if(VERBOSE){cout <<  "\tIntegral: " << Integral;}
   if(Integral>MIN_FIT_COUNTS) {
   
      //-------------------------------------------------------------//
      // First find peaks and identify the first two lines           //
      //-------------------------------------------------------------//
   
      NumPeaks = Spec->Search(Histo, SEARCH_SIGMA, "new",SEARCH_THRESH);
      PeakPositions = Spec->GetPositionX();
      if(VERBOSE){
         cout << "\tPeaks: " << NumPeaks << endl;
         cout << "Commencing peak search for " << En1 << " keV and " << En2 << " keV (Ratio = " << IdealRatio << ")" << endl;
      }
       
      // Loop peaks and identify En1/En2
      PeakFound = 0;
      BestDiff = 1000.0;
      if(NumPeaks>1) {
         for(Peak1=0; Peak1<NumPeaks; Peak1++) {
            if(VERBOSE){cout << "Peak " << Peak1 <<  " - Cent: " << PeakPositions[Peak1] << endl; }
            for(Peak2=0;Peak2<NumPeaks;Peak2++) {
               if(Peak1 != Peak2) {
                  Ratio = PeakPositions[Peak1] / PeakPositions[Peak2];
                  Diff = fabs(Ratio - IdealRatio);
                  //cout << "Calc Gain: " << PEAK_EN2/PeakPos[Peak2] << " Ratio: " << Ratio << " Diff: " << Diff << " Best: " << BestDiff << endl;
                  if(Diff < BestDiff) {  // best match so far                     
                     if( (En2/(PeakPositions[Peak2]/INTEGRATION)) > MIN_GAIN && (En2/(PeakPositions[Peak2]/INTEGRATION)) < MAX_GAIN) { // Gain is sensible
                        BestDiff = Diff;
                        BestPeak1 = Peak1;
                        BestPeak2 = Peak2;
                        PeakFound = 1;
                     }
                  }
               }
            }     
         }
         if(VERBOSE && PeakFound){cout << "Best Peak1: " << BestPeak1 << " BestPeak2: " << BestPeak2 << " Ratio: " << PeakPositions[BestPeak1] / PeakPositions[BestPeak2] << endl;}
      }
      else{
         return -2;
      }   

      if(PeakFound==0) {
         cout << "No matching peaks found!" << endl;
         return -3;
      }      
      

      if(BestDiff < 0.1) {
         
         //-------------------------------------------------------------//
         // Now fit the first two lines and get approx gain             //
         //-------------------------------------------------------------//     
                  
         if(VERBOSE){cout << "Peaks identified, commencing fit..." << endl << endl;}
         
         TF1 *FitRange[NUM_LINES];  // pointers to functions for fitting "gaus" to each line
         
         for(Line=0; Line<2; Line++) { 
            
            if(Line==0) {Peak=BestPeak1;}
            if(Line==1) {Peak=BestPeak2;}
            
            Min = PeakPositions[Peak]-FIT_WIDTH;
            Max = PeakPositions[Peak]+FIT_WIDTH;
            
            if(VERBOSE) {
               cout << "Fitting line " << Line << " Min: " << Min << " Max: " << Max << endl;
            }
            
            if(FIT_BACKGROUND==0) {
               FitRange[Line] =  new TF1("GausFit", "gaus", Min, Max);
               //FitRange[Line]->SetParLimits(2,GAUS_SIGMA_MIN,GAUS_SIGMA_MAX);
            }
            else {
               FitRange[Line] = new TF1 ("GausFlatBack","([0]*exp(-0.5*((x-[1])/[2])**2))+[3]",Min,Max);
               FitRange[Line]->SetParName(0,"Const");
               FitRange[Line]->SetParName(1,"Mean");
               FitRange[Line]->SetParName(2,"Sigma");
               FitRange[Line]->SetParName(3,"Background");
               FitRange[Line]->SetParameter(0,GAUS_CONST_INITIAL);
               FitRange[Line]->SetParameter(1,PeakPositions[Peak]);
               InitialSigma = GAUS_SIGMA_ZERO + ((Sources[Source][Line] / 1000.0)*GAUS_SIGMA_1MEV);
               if(VERBOSE) {cout << "Initial Sigma: " << InitialSigma << endl;}
               FitRange[Line]->SetParameter(2,InitialSigma);
               MinBin = int(Min * CHARGE_BINS / CHARGE_MAX);
               MaxBin = int(Max * CHARGE_BINS / CHARGE_MAX);
               BackBins = int(BACK_WIDTH  * CHARGE_BINS / CHARGE_MAX);
               BackEst = (Histo->Integral(MinBin,(MinBin+BackBins)) + Histo->Integral(MaxBin-BackBins,MaxBin)) / (2*BackBins);
               if(VERBOSE) {cout << "Back Est: " << BackEst << endl;}
               FitRange[Line]->SetParameter(3,BackEst);
               //FitRange[Line]->SetParLimits(2,GAUS_SIGMA_MIN,GAUS_SIGMA_MAX);
            }
            
            if(Line==0) {
               Opts = FitOptions;
            }
            else {
               Opts = FitOptions + "+";   
            }
            
            if(VERBOSE) {cout << "Fit options: " << Opts << endl;}
            Histo->Fit(FitRange[Line],Opts.c_str());  
            //Histo->Fit(FitRange[Line],"R,Q");         
            
            FitRes[Line].Energy = Sources[Source][Line];
            FitRes[Line].Const = FitRange[Line]->GetParameter(0);
            FitRes[Line].dConst = FitRange[Line]->GetParError(0);
            FitRes[Line].Mean = FitRange[Line]->GetParameter(1);
            FitRes[Line].dMean = FitRange[Line]->GetParError(1);
            FitRes[Line].Sigma = FitRange[Line]->GetParameter(2);
            FitRes[Line].dSigma = FitRange[Line]->GetParError(2);
            FitRes[Line].ChiSq = FitRange[Line]->GetChisquare();
            FitRes[Line].NDF = FitRange[Line]->GetNDF();
            
            if(PlotOn || VERBOSE) {
               cout << "Peak " << Line << " Params: " << FitRes[Line].Const << " " << FitRes[Line].Mean << " " << FitRes[Line].Sigma << endl;
               cout << "Peak " << Line << " Errors: " << FitRes[Line].dConst  << " " << FitRes[Line].dMean  << " " << FitRes[Line].dSigma << endl;    
               if(FIT_BACKGROUND==1) {
                  cout << "Background = " << FitRange[Line]->GetParameter(3) << endl;
               }
               cout << "ChiSq: " << FitRes[Line].ChiSq << " NDF: " << FitRes[Line].NDF << " CSPD: " << FitRes[Line].ChiSq/FitRes[Line].NDF << endl;               
            }
            
            if(PlotOn) {            
               c1->cd(1);               
               Histo->Draw();
               c1->Update();
            }
         }
        
         //-------------------------------------------------------------//
         // Calculate a rough calibration                               //
         //-------------------------------------------------------------//
         Chg1 = FitRange[0]->GetParameter(1) / INTEGRATION;
         Chg2 = FitRange[1]->GetParameter(1) / INTEGRATION;
         dChg1 = FitRange[0]->GetParError(1) / INTEGRATION;
         dChg2 = FitRange[0]->GetParError(1) / INTEGRATION;
         
         G = (Sources[Source][1] - Sources[Source][0]) / (Chg2 - Chg1 );
         //*Gain = G;
         dG = sqrt( pow(dChg1,2) + pow(dChg2,2) ) * (1.0/pow(Chg2 - Chg1,2)) * (Sources[Source][1] - Sources[Source][0]);
         //*dGain = dG;
         
         O = Sources[Source][1] - (G * Chg2);
         //*Offset = O;
         dO = sqrt(pow(dG/G,2)+pow(FitRange[1]->GetParError(1)/FitRange[1]->GetParameter(1),2)) * O;
         //*dOffset = dO;
         
         Fit->LinGain[0] = O;
         Fit->dLinGain[0] = dO;
         Fit->LinGain[1] = G;
         Fit->dLinGain[1] = dG;
         
         if(VERBOSE) {
            cout << "Gain: " << G << " +/- " << dG << ", Offset: " << O << " +/- " << dO << endl;
         }
             
         if(PlotOn){ // && Clover==9) {
            c1->cd(1);
            c1->Update();
            App->Run(1);  
         }        
   
         //-------------------------------------------------------------//
         // Loop the remaining lines and fit them                       //
         //-------------------------------------------------------------//
         if(NUM_LINES>2) {
            for(Line=2; Line<NUM_LINES; Line++) {
               
               FitCentre = ((Sources[Source][Line] - O)*INTEGRATION) / G;
               Min = FitCentre - FIT_WIDTH;
               Max = FitCentre + FIT_WIDTH;
               
               if(VERBOSE) {
                  cout << "Fitting line " << Line << " (" << Sources[Source][Line] << " keV), Min: " << Min << " Max: " << Max << endl; 
               }
               
               if(FIT_BACKGROUND==0) {
                  FitRange[Line] =  new TF1("GausFit", "gaus", Min, Max);
                  //FitRange[Line]->SetParLimits(2,GAUS_SIGMA_MIN,GAUS_SIGMA_MAX);
               }
               else {
                  FitRange[Line] = new TF1 ("GausFlatBack","([0]*exp(-0.5*((x-[1])/[2])**2))+[3]",Min,Max);
                  FitRange[Line]->SetParName(0,"Const");
                  FitRange[Line]->SetParName(1,"Mean");
                  FitRange[Line]->SetParName(2,"Sigma");
                  FitRange[Line]->SetParName(3,"Background");
                  FitRange[Line]->SetParameter(0,GAUS_CONST_INITIAL);
                  FitRange[Line]->SetParameter(1,FitCentre);
                  InitialSigma = GAUS_SIGMA_ZERO + ((Sources[Source][Line] / 1000.0)*GAUS_SIGMA_1MEV);
                  if(VERBOSE) {cout << "Initial Sigma: " << InitialSigma << endl;}
                  FitRange[Line]->SetParameter(2,InitialSigma);
                  MinBin = int(Min * CHARGE_BINS / CHARGE_MAX);
                  MaxBin = int(Max * CHARGE_BINS / CHARGE_MAX);
                  BackBins = int(BACK_WIDTH  * CHARGE_BINS / CHARGE_MAX);
                  BackEst = (Histo->Integral(MinBin,(MinBin+BackBins)) + Histo->Integral(MaxBin-BackBins,MaxBin)) / (2*BackBins);
                  if(VERBOSE) {cout << "Background Estimate: " << BackEst << endl;}
                  FitRange[Line]->SetParameter(3,BackEst);
                  //FitRange[Line]->SetParLimits(2,GAUS_SIGMA_MIN,GAUS_SIGMA_MAX);
               }
                            
               //Histo->Fit(FitRange[Line],"R,Q");
               Histo->Fit(FitRange[Line],Opts.c_str());                  
               FitRes[Line].Energy = Sources[Source][Line];
               FitRes[Line].Const = FitRange[Line]->GetParameter(0);
               FitRes[Line].dConst = FitRange[Line]->GetParError(0);
               FitRes[Line].Mean = FitRange[Line]->GetParameter(1);
               FitRes[Line].dMean = FitRange[Line]->GetParError(1);
               FitRes[Line].Sigma = FitRange[Line]->GetParameter(2);
               FitRes[Line].dSigma = FitRange[Line]->GetParError(2);      
               FitRes[Line].ChiSq = FitRange[Line]->GetChisquare();
               FitRes[Line].NDF = FitRange[Line]->GetNDF();
               
               if(PlotOn || VERBOSE) {
                  cout << "Peak " << Line << " Params: " << FitRange[Line]->GetParameter(0) << " " << FitRes[Line].Mean << " " << FitRes[Line].Sigma << endl;
                  cout << "Peak " << Line << " Errors: " << FitRes[Line].dConst  << " " << FitRes[Line].dMean  << " " << FitRes[Line].dSigma << endl;    
                  if(FIT_BACKGROUND==1) {
                     cout << "Background = " << FitRange[Line]->GetParameter(3) << endl;
                  }
                  cout << "ChiSq: " << FitRes[Line].ChiSq << " NDF: " << FitRes[Line].NDF << " CSPD: " << FitRes[Line].ChiSq/FitRes[Line].NDF << endl; 
               }
            }   
                 
            if(PlotOn){ // && Clover==9) {
               
               c1->cd(1);               
               Histo->Draw();
               App->Run(1);  
            } 
         }
         
         //-------------------------------------------------------------//
         // Finally fit energy vs peak centroid for final calibration   //
         //-------------------------------------------------------------//
         if(VERBOSE) {
            cout << endl << "Calibrating..." << endl;
         }
         memset(&Energies,   0.0, sizeof(float)*NUM_LINES);
         memset(&dEnergies,  0.0, sizeof(float)*NUM_LINES);
         memset(&Centroids,  0.0, sizeof(float)*NUM_LINES);
         memset(&dCentroids, 0.0, sizeof(float)*NUM_LINES);
         
         LinesUsed=0;
         for(Line=0; Line<NUM_LINES; Line++) {
            if(FitRes[Line].Const > GAUS_HEIGHT_MIN) {
               if(fabs(FitRes[Line].Sigma) > GAUS_SIGMA_MIN) {
                  if((FitRes[Line].ChiSq / FitRes[Line].NDF) < GAUS_CSPD_MAX) {
                     Energies[LinesUsed] = Sources[Source][Line];
                     Centroids[LinesUsed] = FitRes[Line].Mean / INTEGRATION;
                     dCentroids[LinesUsed] = FitRes[Line].dMean / INTEGRATION;
                     if(VERBOSE) {
                        cout << Energies[LinesUsed] << " keV\t" << Centroids[LinesUsed] << " +/- " << dCentroids[LinesUsed] << " ch" << endl;
                     }
                     // Pass fit results back out to main
                     Fit->PeakFits[Line].Energy = FitRes[Line].Energy;    
                     Fit->PeakFits[Line].Const = FitRes[Line].Const ;
                     Fit->PeakFits[Line].dConst = FitRes[Line].dConst;
                     Fit->PeakFits[Line].Mean = FitRes[Line].Mean   ;
                     Fit->PeakFits[Line].dMean = FitRes[Line].dMean  ;
                     Fit->PeakFits[Line].Sigma = FitRes[Line].Sigma  ;
                     Fit->PeakFits[Line].dSigma = FitRes[Line].dSigma  ;     
                     Fit->PeakFits[Line].ChiSq = FitRes[Line].ChiSq  ;
                     Fit->PeakFits[Line].NDF = FitRes[Line].NDF     ;  
                     Fit->FitSuccess[Line]=1;              
                     LinesUsed += 1;
                  }
               }   
            }         
         }
         
         TGraphErrors CalibPlot(LinesUsed,Centroids,Energies,dCentroids,dEnergies);
         if(VERBOSE) {
            cout << LinesUsed << " lines used for calibration ";
            for(i=0; i<LinesUsed; i++) {
               cout << Energies[i] << " ";
            }
            cout << endl;
         }
         
         TF1 *CalibFitLin  = new TF1("CalibFitLin","[0] + ([1]*x)",0.0,CHARGE_MAX);
         CalibFitLin->SetParName(0,"Offset");
         CalibFitLin->SetParName(1,"Gain");
         CalibFitLin->SetParameter(0,0.0);
         CalibFitLin->SetParameter(1,INITIAL_GAIN);
         CalibFitLin->SetLineColor(4); // Blue?
         TF1 *CalibFitQuad = new TF1 ("CalibFitQuad","[0] + ([1]*x) + ([2]*x*x)",0.0,CHARGE_MAX);
         CalibFitQuad->SetParName(0,"Offset");
         CalibFitQuad->SetParName(1,"Gain");
         CalibFitQuad->SetParName(2,"Quad");
         CalibFitQuad->SetParameter(0,0.0);
         CalibFitQuad->SetParameter(1,INITIAL_GAIN);
         CalibFitQuad->SetParameter(2,0.0);
         CalibFitQuad->SetLineColor(2); //Red?
         
         Opts = FitOptions;
         if(LinesUsed > 1) {        
            CalibPlot.Fit(CalibFitLin,Opts.c_str());   
            Fit->LinGainFit[0] = CalibFitLin->GetParameter(0);
            Fit->dLinGainFit[0] = CalibFitLin->GetParError(0);
            Fit->LinGainFit[1] = CalibFitLin->GetParameter(1);
            Fit->dLinGainFit[1] = CalibFitLin->GetParError(1);
            Fit->LinGainFit[2] = CalibFitLin->GetChisquare() / CalibFitLin->GetNDF();
         }
         Opts += "+";   
         if(LinesUsed > 2) {
            CalibPlot.Fit(CalibFitQuad,Opts.c_str());
            Fit->QuadGainFit[0] = CalibFitQuad->GetParameter(0);
            Fit->dQuadGainFit[0] = CalibFitQuad->GetParError(0);
            Fit->QuadGainFit[1] = CalibFitQuad->GetParameter(1);
            Fit->dQuadGainFit[1] = CalibFitQuad->GetParError(1);
            Fit->QuadGainFit[2] = CalibFitQuad->GetParameter(2);
            Fit->dQuadGainFit[2] = CalibFitQuad->GetParError(2);
            Fit->QuadGainFit[3] = CalibFitQuad->GetChisquare() / CalibFitQuad->GetNDF();
         }
         
         if(VERBOSE) {
            cout << endl << "Two point solution:" << endl;
            cout <<  "Offset = " << O << " +/- " << dO << "\tGain = " << G << " +/- " << dG << endl << endl;
            cout << "Linear fit: " << endl;
            cout << "Offset = " << CalibFitLin->GetParameter(0) << " +/- " <<  CalibFitLin->GetParError(0) << "\t";         
            cout << "Gain = " << CalibFitLin->GetParameter(1) << " +/- " <<  CalibFitLin->GetParError(1) << endl;
            cout << "CSPD = " << CalibFitLin->GetChisquare() / CalibFitLin->GetNDF() << endl << endl;
            cout << "Quadratic fit: " << endl;
            cout << "Offset = " << CalibFitQuad->GetParameter(0) << " +/- " <<  CalibFitQuad->GetParError(0) << "\t";            
            cout << "Gain = " << CalibFitQuad->GetParameter(1) << " +/- " <<  CalibFitQuad->GetParError(1) << "\t";
            cout << "Quad = " << CalibFitQuad->GetParameter(2) << " +/- " <<  CalibFitQuad->GetParError(2) << endl;
            cout << "CSPD = " << CalibFitQuad->GetChisquare() / CalibFitQuad->GetNDF() << endl << endl;
         }
         
         
         
         if(PLOT_CALIB && PlotOn) {
            c1->cd(2);
            CalibPlot.SetMarkerColor(2); 
            CalibPlot.SetMarkerStyle(20); 
            CalibPlot.SetMarkerSize(1.0);
            CalibPlot.Draw("AP");           
            App->Run(1);  
         }
         
      }
      else {
         if(VERBOSE){
            cout << "Ratio of peak centroids (" << BestDiff << ") does not match the expected value (" << IdealRatio << ")" << endl;
         }
         return -4;
      }   
   }
   
   else {
      if(VERBOSE){
         cout << "\tNot enough counts to fit (" << Integral << ")" << endl;        
      }
      return -1;
   }
   
   return 1;

}



