// Fit spectrum function for use with Calib and CalibOffline codes

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
//#include "TTigFragment.h"
//#include "TFSPC_Info.h"
//#include "TSharc.h"
//#include "TTigress.h"
//#include "TRf.h"
//#include "TTriFoil.h"

// My libraries
#include "Calib.h"
#include "Main.h"
//#include "FitGammaSpectrum.h"

extern TApplication* App;

int FitGammaSpectrum(TH1F* Histo, SpectrumFit *Fit, int Source, int Integration, int PlotOn ) {

   int Peak, Peak1, Peak2, NumPeaks, BestPeak1, BestPeak2, PeakFound;  // for looping peaks and finding correct ones
   float Ratio, Diff, BestDiff;  // test quality of peak match
   float Min, Max, BackEst;
   int i;
   int MinBin, MaxBin, BackBins;
   Float_t* PeakPositions;  // for output of root peak search
   int Line, LinesUsed;   
   float FitCentre, G, O, dG, dO;  // Fits and values for fast 2 point calibration.
   float Energies[NUM_LINES+1], dEnergies[NUM_LINES+1], Centroids[NUM_LINES+1], dCentroids[NUM_LINES+1]; // Data for full calibration
   float En1, En2;   // line energies for two point calib    
   En1 = Sources[Source][0];
   En2 = Sources[Source][1];
   float IdealRatio = En1/En2;  
   float Chg1, Chg2, dChg1, dChg2;  // charge and erros for 2 point         
   FitResult FitRes[NUM_LINES];  // store full fit results
   int Integral;  // counts in spectrum
   TSpectrum *Spec = new TSpectrum();
   
   //float InitialSigma = 0.0;
   float InitialSigma = 50.0;
   float InitialGain = 0.16;
   
   // Fitting stuff
   std::string FitOptions = ("RQEM");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
   std::string Opts;
   
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
                  cout << "Calc Gain: " << En2/PeakPositions[Peak2] << " Ratio: " << Ratio << " Diff: " << Diff << " Best: " << BestDiff << endl;
                  if(Diff < BestDiff) {  // best match so far                     
                     //if( (En2/(PeakPositions[Peak2]/Integration)) > MIN_GAIN && (En2/(PeakPositions[Peak2]/Integration)) < MAX_GAIN) { // Gain is sensible
                        BestDiff = Diff;
                        BestPeak1 = Peak1;
                        BestPeak2 = Peak2;
                        PeakFound = 1;
                     //}
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
            
            Min = float(PeakPositions[Peak]-FIT_WIDTH);
            Max = float(PeakPositions[Peak]+FIT_WIDTH);
            
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
               cCalib1->cd(1);               
               Histo->Draw();
               cCalib1->Update();
            }
         }
        
         //-------------------------------------------------------------//
         // Calculate a rough calibration                               //
         //-------------------------------------------------------------//
         Chg1 = FitRange[0]->GetParameter(1) / Integration;
         Chg2 = FitRange[1]->GetParameter(1) / Integration;
         dChg1 = FitRange[0]->GetParError(1) / Integration;
         dChg2 = FitRange[0]->GetParError(1) / Integration;
         
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
            cCalib1->cd(1);
            cCalib1->Update();
            App->Run(1);  
         }        
   
         //-------------------------------------------------------------//
         // Loop the remaining lines and fit them                       //
         //-------------------------------------------------------------//
         if(NUM_LINES>2) {
            for(Line=2; Line<NUM_LINES; Line++) {
               
               FitCentre = ((Sources[Source][Line] - O)*Integration) / G;
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
               
               cCalib1->cd(1);               
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
         memset(&Energies,   0.0, sizeof(float)*(NUM_LINES+1));
         memset(&dEnergies,  0.0, sizeof(float)*(NUM_LINES+1));
         memset(&Centroids,  0.0, sizeof(float)*(NUM_LINES+1));
         memset(&dCentroids, 0.0, sizeof(float)*(NUM_LINES+1));
         
         LinesUsed=0;
         for(Line=0; Line<NUM_LINES; Line++) {
            if(FitRes[Line].Const > GAUS_HEIGHT_MIN) {
               if(fabs(FitRes[Line].Sigma) > GAUS_SIGMA_MIN) {
                  if((FitRes[Line].ChiSq / FitRes[Line].NDF) < GAUS_CSPD_MAX) {
                     Energies[LinesUsed] = Sources[Source][Line];
                     Centroids[LinesUsed] = FitRes[Line].Mean / Integration;
                     dCentroids[LinesUsed] = FitRes[Line].dMean / Integration;
                     if(VERBOSE) {
                        cout << Energies[LinesUsed] << " keV\t" << Centroids[LinesUsed] << " +/- " << dCentroids[LinesUsed] << " ch" << endl;
                     }
                     // Pass fit results back out to main
                     Fit->PeakFits[Line].Energy = FitRes[Line].Energy;    
                     Fit->PeakFits[Line].Const = FitRes[Line].Const;
                     Fit->PeakFits[Line].dConst = FitRes[Line].dConst;
                     Fit->PeakFits[Line].Mean = FitRes[Line].Mean;
                     Fit->PeakFits[Line].dMean = FitRes[Line].dMean;
                     Fit->PeakFits[Line].Sigma = FitRes[Line].Sigma;
                     Fit->PeakFits[Line].dSigma = FitRes[Line].dSigma;     
                     Fit->PeakFits[Line].ChiSq = FitRes[Line].ChiSq;
                     Fit->PeakFits[Line].NDF = FitRes[Line].NDF;  
                     Fit->FitSuccess[Line]=1;              
                     LinesUsed += 1;
                  }
               }   
            }         
         }
         if(INCLUDE_ZERO) {
            Energies[LinesUsed] = 0.0;
            Centroids[LinesUsed] = 0.0;
            dCentroids[LinesUsed] = ZERO_ERR;
            if(VERBOSE) {
               cout << Energies[LinesUsed] << " keV\t" << Centroids[LinesUsed] << " +/- " << dCentroids[LinesUsed] << " ch" << endl;
            }
            Fit->PeakFits[NUM_LINES].Energy = 0.0;    
            Fit->PeakFits[NUM_LINES].Const = 0.0;
            Fit->PeakFits[NUM_LINES].dConst = 0.0;
            Fit->PeakFits[NUM_LINES].Mean = 0.0;
            Fit->PeakFits[NUM_LINES].dMean = 1.0;
            Fit->PeakFits[NUM_LINES].Sigma = 0.0;
            Fit->PeakFits[NUM_LINES].dSigma = 0.0;     
            Fit->PeakFits[NUM_LINES].ChiSq = 0.0;
            Fit->PeakFits[NUM_LINES].NDF = 1;  
            Fit->FitSuccess[Line]=1;              
            LinesUsed += 1;
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
            cCalib1->cd(2);
            CalibPlot.SetMarkerColor(2); 
            CalibPlot.SetMarkerStyle(20); 
            CalibPlot.SetMarkerSize(1.0);
            CalibPlot.SetTitle("Calibration");
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
