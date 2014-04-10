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
#include "CalibTools.h"
//#include "FitGammaSpectrum.h"

extern TApplication *App;

int FitGammaSpectrum(TH1F * Histo, SpectrumFit * Fit, FitSettings Settings)
{

   // Histo:         pointer to histo to be fitted
   // Fit:           SpectrumFit structure into which the fit results will be written
   // Settings:
   //    Source:        Number of source to be fitted
   //    Integration:   Integration used in charge evaluation (set to 1 if no integration or if this is the waveform derived charge)
   //    Dispersion:    number of bins in charge spectrum per unit
   //    PlotOn:        should these fits be plotted

   // Variables needed:
   int Peak, Peak1, Peak2, NumPeaks, BestPeak1, BestPeak2, PeakFound;   // for looping peaks and finding correct ones
   int i, temp;
   float Ratio, Diff, BestDiff; // test quality of peak match
   float Centre;
   Float_t *PeakPositions;      // for output of root peak search
   int Line, LinesUsed;
   float FitCentre, G, O, dG, dO;       // Fits and values for fast 2 point calibration.
   float Energies[MAX_LINES + 1], dEnergies[MAX_LINES + 1], Centroids[MAX_LINES + 1], dCentroids[MAX_LINES + 1];        // Data for full calibration
   float En1, En2;              // line energies for two point calib
   En1 = Config.Sources[Settings.Source][0];
   En2 = Config.Sources[Settings.Source][1];
   float IdealRatio = En1 / En2;
   float Chg1, Chg2, dChg1, dChg2;      // charge and erros for 2 point
   FitResult FitRes[MAX_LINES]; // store full fit results
   int Integral;                // counts in spectrum
   TSpectrum *Spec = new TSpectrum();
   TF1 *FitRange[MAX_LINES];    // pointers to functions for fitting each line
   float InitialGain = 0.16;
   // Fitting stuff
   std::string FitOptions = ("RQEM");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
   std::string Opts;



   if (Histo) {
      Integral = Histo->Integral();
   } else {
      cout << "Bad histo! " << endl;
   }
   if (Config.PrintVerbose) {
      cout << "\tIntegral: " << Integral;
   }
   if (Integral > MIN_FIT_COUNTS) {

      if (PLOT_FITS) {
         cCalib1->cd(1);
         cCalib1->Modified();
         cCalib1->Update();
      }
      //-------------------------------------------------------------//
      // First find peaks and identify the first two lines           //
      //-------------------------------------------------------------//

      NumPeaks = Spec->Search(Histo, Settings.SearchSigma, "new", Settings.SearchThresh);
      PeakPositions = Spec->GetPositionX();
      if (Config.PrintVerbose) {
         cout << "\tPeaks: " << NumPeaks << endl;
         cout << "Commencing peak search for " << En1 << " keV and " << En2 << " keV (Ratio = " << IdealRatio << ")" <<
             endl;
      }
      // Loop peaks and identify En1/En2
      PeakFound = 0;
      BestDiff = 1000.0;
      if (NumPeaks > 1) {
         for (Peak1 = 0; Peak1 < NumPeaks; Peak1++) {
            if (Config.PrintVerbose) {
               cout << "Peak " << Peak1 << " - Cent: " << PeakPositions[Peak1] << endl;
            }
            for (Peak2 = 0; Peak2 < NumPeaks; Peak2++) {
               if (Peak1 != Peak2) {
                  Ratio = PeakPositions[Peak1] / PeakPositions[Peak2];
                  Diff = fabs(Ratio - IdealRatio);
                  //cout << "Calc Gain: " << En2 /
                  //  PeakPositions[Peak2] << " Ratio: " << Ratio << " Diff: " << Diff << " Best: " << BestDiff << endl;
                  if (Diff < BestDiff) {        // best match so far
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
         if (Config.PrintVerbose && PeakFound) {
            cout << "Best Peak1: " << BestPeak1 << " BestPeak2: " << BestPeak2 << " Ratio: " << PeakPositions[BestPeak1]
                / PeakPositions[BestPeak2] << endl;
         }
      } else {
         return -2;
      }

      if (PeakFound == 0) {
         cout << "No matching peaks found!" << endl;
         return -3;
      }


      if (BestDiff < 0.1) {

         //-------------------------------------------------------------//
         // Now fit the first two lines and get approx gain             //
         //-------------------------------------------------------------//

         if (Config.PrintVerbose) {
            cout << "Peaks identified, commencing fit..." << endl << endl;
         }

         for (Line = 0; Line < 2; Line++) {
            if (Line == 0) {
               Peak = BestPeak1;
            }
            if (Line == 1) {
               Peak = BestPeak2;
            }

            Centre = PeakPositions[Peak];
            if (Config.PrintVerbose) {
               cout << "Line: " << Line << " (Energy = " << Config.Sources[Settings.Source][Line] << " keV)" << endl;
               cout << "-------------------------------------------------------------" << endl << endl;
            }
            FitSinglePeak(Histo, Line, Centre, FitRange[Line], &FitRes[Line], Settings);

         }

         //cin >> temp;

         //-------------------------------------------------------------//
         // Calculate a rough calibration                               //
         //-------------------------------------------------------------//

         /*cout << "Here! " << endl;
            cout << "Set.Int " << Settings.Integration << endl;
            cout << "FitRes[0].Mean " << FitRes[0].Mean << endl;
            cout << "FitRes[1].Mean " << FitRes[1].Mean << endl; */

         Chg1 = FitRes[0].Mean / Settings.Integration;
         Chg2 = FitRes[1].Mean / Settings.Integration;
         dChg1 = FitRes[0].dMean / Settings.Integration;
         dChg2 = FitRes[1].dMean / Settings.Integration;

         /*cout << "Chg1 " << Chg1 << endl;
            cout << "Chg2 " << Chg2 << endl;

            cout << "Source1 " << Sources[Settings.Source][0] << endl;
            cout << "Source2 " << Sources[Settings.Source][1] << endl; */

         // Gain = (E2 - E1) / (CHG2 - CHG1)       
         G = (Config.Sources[Settings.Source][1] - Config.Sources[Settings.Source][0]) / (Chg2 - Chg1);
         // dGain = sqrt( dCHG1^2 + cCHG2^2 )  *  (E2 - E1)  /  (CHG2 - CHG1)^2     [assuming no error on E]
         dG = sqrt(pow(dChg1, 2) + pow(dChg2, 2)) * (1.0 / pow(Chg2 - Chg1, 2)) * (Config.Sources[Settings.Source][1] -
                                                                                   Config.Sources[Settings.Source][0]);
         // Offset = E2 - (Gain * CHG2)
         O = Config.Sources[Settings.Source][1] - (G * Chg2);
         // dOffset = Offset * sqrt( (dCHG2/CHG2)^2 + (dGain/Gain)^2 )
         dO = sqrt(pow(dG / G, 2) + pow(dChg2 / Chg2, 2)) * O;
         //*dOffset = dO;

         Fit->LinGain[0] = O;
         Fit->dLinGain[0] = dO;
         Fit->LinGain[1] = G;
         Fit->dLinGain[1] = dG;

         if (Config.PrintVerbose) {
            cout << "Gain: " << G << " +/- " << dG << ", Offset: " << O << " +/- " << dO << endl;
         }

         if (Settings.PlotOn) { // && Clover==9) {
            cCalib1->cd(1);
            cCalib1->Update();
            App->Run(1);
            //App->Run();
         }
         //-------------------------------------------------------------//
         // Loop the remaining lines and fit them                       //
         //-------------------------------------------------------------//
         if (Config.Sources[Settings.Source].size() > 2) {
            for (Line = 2; Line < Config.Sources[Settings.Source].size(); Line++) {

               if (Config.PrintVerbose) {
                  cout << "Line: " << Line << " (Energy = " << Config.Sources[Settings.Source][Line] << " keV)" << endl;
                  cout << "-------------------------------------------------------------" << endl << endl;
               }
               Centre = ((Config.Sources[Settings.Source][Line] - O) / G) * Settings.Integration;       // Get centre from energy and initial calibration
               FitSinglePeak(Histo, Line, Centre, FitRange[Line], &FitRes[Line], Settings);
            }
         }
         //cin >> temp;

         //-------------------------------------------------------------//
         // Finally fit energy vs peak centroid for final calibration   //
         //-------------------------------------------------------------//
         if (Config.PrintVerbose) {
            cout << endl << "Calibrating..." << endl;
         }
         if (PLOT_FITS || PLOT_CALIB || PLOT_CALIB_SUMMARY || PLOT_RESIDUAL) {
            cCalib1a->cd(1);
         }

         memset(&Energies, 0.0, sizeof(float) * (MAX_LINES + 1));
         memset(&dEnergies, 0.0, sizeof(float) * (MAX_LINES + 1));
         memset(&Centroids, 0.0, sizeof(float) * (MAX_LINES + 1));
         memset(&dCentroids, 0.0, sizeof(float) * (MAX_LINES + 1));

         LinesUsed = 0;
         for (Line = 0; Line < Config.Sources[Settings.Source].size(); Line++) {

            //  Check on these conditions:

            //cout << "Here!" << endl;

            if (FitRes[Line].Const > GAUS_HEIGHT_MIN) {
               //cout << "Here!!" << endl;
               if ((FitRes[Line].ChiSq / FitRes[Line].NDF) < GAUS_CSPD_MAX) {
                  //cout << "Here!!!" << endl;
                  Energies[LinesUsed] = Config.Sources[Settings.Source][Line];
                  Centroids[LinesUsed] = FitRes[Line].Mean / Settings.Integration;
                  dCentroids[LinesUsed] = FitRes[Line].dMean / Settings.Integration;
                  if (Config.PrintVerbose) {
                     cout << Energies[LinesUsed] << " keV\t" << Centroids[LinesUsed] << " +/- " <<
                         dCentroids[LinesUsed] << " ch" << endl;
                  }
                  // Pass fit results back out to main                     
                  memcpy(&Fit->PeakFits[Line], &FitRes[Line], sizeof(FitResult));
                  Fit->FitSuccess[Line] = 1;

                  /*cout << "Mean: " << Fit->PeakFits[Line].Mean << " " << FitRes[Line].Mean << endl;
                     cout << "Sigma: " << Fit->PeakFits[Line].Sigma << " " << FitRes[Line].Sigma << endl;
                     cout << "NDF: " << Fit->PeakFits[Line].NDF << " " << FitRes[Line].NDF << endl; */

                  LinesUsed += 1;
               }
               //}
            }
         }
         if (Settings.FitZero) {
            Energies[LinesUsed] = 0.0;
            Centroids[LinesUsed] = 0.0;
            dCentroids[LinesUsed] = ZERO_ERR;
            if (Config.PrintVerbose) {
               cout << Energies[LinesUsed] << " keV\t" << Centroids[LinesUsed] << " +/- " << dCentroids[LinesUsed] <<
                   " ch" << endl;
            }
            Fit->PeakFits[Config.Sources[Settings.Source].size()].Energy = 0.0;
            Fit->PeakFits[Config.Sources[Settings.Source].size()].Const = 0.0;
            Fit->PeakFits[Config.Sources[Settings.Source].size()].dConst = 0.0;
            Fit->PeakFits[Config.Sources[Settings.Source].size()].Mean = 0.0;
            Fit->PeakFits[Config.Sources[Settings.Source].size()].dMean = 1.0;
            Fit->PeakFits[Config.Sources[Settings.Source].size()].Sigma = 0.0;
            Fit->PeakFits[Config.Sources[Settings.Source].size()].dSigma = 0.0;
            Fit->PeakFits[Config.Sources[Settings.Source].size()].ChiSq = 0.0;
            Fit->PeakFits[Config.Sources[Settings.Source].size()].NDF = 1;
            Fit->FitSuccess[Line + 1] = 1;
            LinesUsed += 1;
         }
         Fit->LinesUsed = LinesUsed;
         TGraphErrors CalibPlot(LinesUsed, Centroids, Energies, dCentroids, dEnergies);
         if (Config.PrintVerbose) {
            cout << LinesUsed << " lines used for calibration ";
            for (i = 0; i < LinesUsed; i++) {
               cout << Energies[i] << " ";
            }
            cout << endl;
         }

         TF1 *CalibFitLin = new TF1("CalibFitLin", "[0] + ([1]*x)", 0.0, CHARGE_MAX);
         CalibFitLin->SetParName(0, "Offset");
         CalibFitLin->SetParName(1, "Gain");
         CalibFitLin->SetParameter(0, 0.0);
         CalibFitLin->SetParameter(1, INITIAL_GAIN);
         CalibFitLin->SetLineColor(4);  // Blue?
         TF1 *CalibFitQuad = new TF1("CalibFitQuad", "[0] + ([1]*x) + ([2]*x*x)", 0.0, CHARGE_MAX);
         CalibFitQuad->SetParName(0, "Offset");
         CalibFitQuad->SetParName(1, "Gain");
         CalibFitQuad->SetParName(2, "Quad");
         CalibFitQuad->SetParameter(0, 0.0);
         CalibFitQuad->SetParameter(1, INITIAL_GAIN);
         CalibFitQuad->SetParameter(2, 0.0);
         CalibFitQuad->SetLineColor(2); //Red?

         Opts = FitOptions;
         if (LinesUsed > 1) {
            CalibPlot.Fit(CalibFitLin, Opts.c_str());
            Fit->LinGainFit[0] = CalibFitLin->GetParameter(0);
            Fit->dLinGainFit[0] = CalibFitLin->GetParError(0);
            Fit->LinGainFit[1] = CalibFitLin->GetParameter(1);
            Fit->dLinGainFit[1] = CalibFitLin->GetParError(1);
            Fit->LinGainFit[2] = CalibFitLin->GetChisquare() / CalibFitLin->GetNDF();
         }
         Opts += "+";
         if (LinesUsed > 2) {
            CalibPlot.Fit(CalibFitQuad, Opts.c_str());
            Fit->QuadGainFit[0] = CalibFitQuad->GetParameter(0);
            Fit->dQuadGainFit[0] = CalibFitQuad->GetParError(0);
            Fit->QuadGainFit[1] = CalibFitQuad->GetParameter(1);
            Fit->dQuadGainFit[1] = CalibFitQuad->GetParError(1);
            Fit->QuadGainFit[2] = CalibFitQuad->GetParameter(2);
            Fit->dQuadGainFit[2] = CalibFitQuad->GetParError(2);
            Fit->QuadGainFit[3] = CalibFitQuad->GetChisquare() / CalibFitQuad->GetNDF();
         }

         if (Config.PrintVerbose) {
            cout << endl << "Two point solution:" << endl;
            cout << "Offset = " << O << " +/- " << dO << "\tGain = " << G << " +/- " << dG << endl << endl;
            cout << "Linear fit: " << endl;
            cout << "Offset = " << CalibFitLin->GetParameter(0) << " +/- " << CalibFitLin->GetParError(0) << "\t";
            cout << "Gain = " << CalibFitLin->GetParameter(1) << " +/- " << CalibFitLin->GetParError(1) << endl;
            cout << "CSPD = " << CalibFitLin->GetChisquare() / CalibFitLin->GetNDF() << endl << endl;
            cout << "Quadratic fit: " << endl;
            cout << "Offset = " << CalibFitQuad->GetParameter(0) << " +/- " << CalibFitQuad->GetParError(0) << "\t";
            cout << "Gain = " << CalibFitQuad->GetParameter(1) << " +/- " << CalibFitQuad->GetParError(1) << "\t";
            cout << "Quad = " << CalibFitQuad->GetParameter(2) << " +/- " << CalibFitQuad->GetParError(2) << endl;
            cout << "CSPD = " << CalibFitQuad->GetChisquare() / CalibFitQuad->GetNDF() << endl << endl;
         }



         if (PLOT_CALIB && Settings.PlotOn) {
            cCalib1a->cd(1);
            CalibPlot.SetMarkerColor(2);
            CalibPlot.SetMarkerStyle(20);
            CalibPlot.SetMarkerSize(1.0);
            CalibPlot.SetTitle("Calibration");
            CalibPlot.Draw("AP");
            App->Run(1);
            //App->Run();
         }

      } else {
         if (Config.PrintVerbose) {
            cout << "Ratio of peak centroids (" << BestDiff << ") does not match the expected value (" << IdealRatio <<
                ")" << endl;
         }
         return -4;
      }
   }

   else {
      if (Config.PrintVerbose) {
         cout << "\tNot enough counts to fit (" << Integral << ")" << endl;
      }
      return -1;
   }

   return LinesUsed;

}


int FitSinglePeak(TH1F * Histo, int Line, float Centre, TF1 * FitRange, FitResult * FitRes, FitSettings Settings)
{

   int MinBin, MaxBin, BackBins, ConstEst;
   float MinkeV, MaxkeV, Min, Max, BackEst;
   int i;
   float SigmaZero = 0.0;
   float Sigma1MeV = 0.0;
   float InitialSigma = 0.0;
   // Fitting stuff
   std::string FitOptions = ("RQEM");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
   std::string Opts;

   // Find min, max, and other things needed for the fit
   // In keV....
   MinkeV = float (Config.Sources[Settings.Source][Line] - FIT_WIDTH_KEV);
   MaxkeV = float (Config.Sources[Settings.Source][Line] + FIT_WIDTH_KEV);
   // In charge/ADC units
   Min = MinkeV * (Centre / Config.Sources[Settings.Source][Line]);
   Max = MaxkeV * (Centre / Config.Sources[Settings.Source][Line]);
   // In bins
   MinBin = int (Min * Settings.Dispersion);
   if (MinBin < 0) {
      MinBin = 0;
   }
   MaxBin = int (Max * Settings.Dispersion);
   // Estimate Const from maximum bin content
   ConstEst = 0;
   for (i = MinBin; i <= MaxBin; i++) {
      if (Histo->GetBinContent(i) > ConstEst) {
         ConstEst = Histo->GetBinContent(i);
      }
   }
   // find number of bins to be used for initial background estimate
   BackBins = int (BACK_WIDTH_KEV * (Centre / Config.Sources[Settings.Source][Line]) * Settings.Dispersion);
   // Find the estimate for sigma
   SigmaZero = (Settings.SigmaEstZero * (Centre / Config.Sources[Settings.Source][Line]));
   Sigma1MeV = (Settings.SigmaEst1MeV * (Centre / Config.Sources[Settings.Source][Line]));
   InitialSigma = SigmaZero + ((Config.Sources[Settings.Source][Line] / 1000.0) * Sigma1MeV);

   /*cout << "MinkeV: " << MinkeV << " Min: " << Min << " MinBin: " << MinBin << endl;
      cout << "MaxkeV: " << MaxkeV << " Max: " << Max << " MaxBin: " << MaxBin << endl;
      cout << "BackBins: " << BackBins << endl;
      cout << "InitialSigma: " << InitialSigma << endl;
      cout << "Initial Const: " << ConstEst << endl;
      cout << "Dispersion: " << Settings.Dispersion << endl;
    */

   if (Config.PrintVerbose) {
      cout << "Fitting line " << Line << " Min: " << Min << " Max: " << Max << endl;
   }
   if (PLOT_FITS) {
      cCalib1->cd();
   }

   if (FIT_BACKGROUND == 0) {
      FitRange = new TF1("GausFit", "gaus", Min, Max);
      //FitRange->SetParLimits(2,GAUS_SIGMA_MIN,GAUS_SIGMA_MAX);
   } else {
      FitRange = new TF1("GausFlatBack", "([0]*exp(-0.5*((x-[1])/[2])**2))+[3]", Min, Max);
      FitRange->SetParName(0, "Const");
      FitRange->SetParName(1, "Mean");
      FitRange->SetParName(2, "Sigma");
      FitRange->SetParName(3, "Constant Background");
      FitRange->SetParameter(0, ConstEst);
      FitRange->SetParameter(1, Centre);
      if (Config.PrintVerbose) {
         cout << "Initial Sigma: " << InitialSigma << endl;
      }
      FitRange->SetParameter(2, InitialSigma);
      BackEst =
          (Histo->Integral(MinBin, (MinBin + BackBins)) + Histo->Integral(MaxBin - BackBins, MaxBin)) / (2 * BackBins);
      if (Config.PrintVerbose) {
         cout << "Back Est: " << BackEst << endl;
      }
      FitRange->SetParameter(3, BackEst);
      //FitRange->SetParLimits(2,GAUS_SIGMA_MIN,GAUS_SIGMA_MAX);
   }

   if (Line == 0) {
      Opts = FitOptions;
   } else {
      Opts = FitOptions + "+";
   }

   if (Config.PrintVerbose) {
      cout << "Fit options: " << Opts << endl;
   }

   Histo->Fit(FitRange, Opts.c_str());
   //Histo->Fit(FitRange,"R,Q");

   FitRes->Energy = Config.Sources[Settings.Source][Line];
   FitRes->Const = FitRange->GetParameter(0);
   FitRes->dConst = FitRange->GetParError(0);
   FitRes->Mean = FitRange->GetParameter(1);
   FitRes->dMean = FitRange->GetParError(1);
   FitRes->Sigma = FitRange->GetParameter(2);
   FitRes->dSigma = FitRange->GetParError(2);
   FitRes->ConstantBG = FitRange->GetParameter(3);
   FitRes->dConstantBG = FitRange->GetParError(3);
   FitRes->ChiSq = FitRange->GetChisquare();
   FitRes->NDF = FitRange->GetNDF();

   if (Settings.PlotOn || Config.PrintVerbose) {
      cout << "Peak " << Line << " Params: " << FitRes->Const << " " << FitRes->Mean << " " << FitRes->Sigma << endl;
      cout << "Peak " << Line << " Errors: " << FitRes->dConst << " " << FitRes->dMean << " " << FitRes->dSigma << endl;
      if (FIT_BACKGROUND == 1) {
         cout << "Background = " << FitRange->GetParameter(3) << endl;
      }
      cout << "ChiSq: " << FitRes->ChiSq << " NDF: " << FitRes->NDF << " CSPD: " << FitRes->ChiSq /
          FitRes->NDF << endl << endl;
   }

   if (Settings.PlotOn) {
      cCalib1->cd();
      Histo->Draw();
      cCalib1->Update();
      //App->Run(1);
      //App->Run();
   }

   return 1;
}


int CalibrationReport(SpectrumFit * Fit, ofstream & ReportOut, std::string HistName, FitSettings Settings)
{

   int i, NumFits;
   float Energies[MAX_LINES], Residuals[MAX_LINES];
   float CalibEn;

   //char HistName[1024];
   //sprintf(HistName,"TempHist");
   
   if(Config.CalReport) {
      // Heading for this channel
      ReportOut << endl << "------------------------------------------" << endl << HistName << endl <<
          "------------------------------------------" << endl << endl;
      // Peak fit info
      ReportOut << Fit->LinesUsed << " lines used in calibration." << endl;
      ReportOut << "Individual fit results:" << endl;
      ReportOut << "(Const +/- err\tMean +/- err\tSigma +/- err\nChiSq\tNDF\tCSPD)" << endl << endl;
      for (i = 0; i <= Config.Sources[Settings.Source].size(); i++) {
         if (Fit->FitSuccess[i] == 1) {
            ReportOut << Fit->PeakFits[i].Energy << " keV" << endl;
            ReportOut << Fit->PeakFits[i].Const << " +/- " << Fit->PeakFits[i].dConst << "\t";
            ReportOut << Fit->PeakFits[i].Mean << " +/- " << Fit->PeakFits[i].dMean << "\t";
            ReportOut << Fit->PeakFits[i].Const << " +/- " << Fit->PeakFits[i].dConst << endl;
            ReportOut << Fit->PeakFits[i].ChiSq << "\t" << Fit->PeakFits[i].NDF << "\t" << Fit->PeakFits[i].ChiSq /
                Fit->PeakFits[i].NDF;
            ReportOut << endl << endl;
         }
      }
      // Calibration...
      ReportOut << "Linear Solution: Offset = " << Fit->LinGain[0] << " +/- " << Fit->dLinGain[0] << "\t";
      ReportOut << "Gain = " << Fit->LinGain[1] << " +/- " << Fit->dLinGain[1] << endl;
      ReportOut << "Linear Fit: Offset = " << Fit->LinGainFit[0] << " +/- " << Fit->dLinGainFit[0] << "\t";
      ReportOut << "Gain = " << Fit->LinGainFit[1] << " +/- " << Fit->dLinGainFit[1] << "\t";
      ReportOut << "CSPD = " << Fit->LinGainFit[2] << endl;
      ReportOut << "Quadratic Fit: Offset = " << Fit->QuadGainFit[0] << " +/- " << Fit->dQuadGainFit[0] << "\t";
      ReportOut << "Gain = " << Fit->QuadGainFit[1] << " +/- " << Fit->dQuadGainFit[1] << "\t";
      ReportOut << "Quad = " << Fit->QuadGainFit[2] << " +/- " << Fit->dQuadGainFit[2] << "\t";
      ReportOut << "CSPD = " << Fit->LinGainFit[3] << endl;
      // Residual from quadratic fit
      ReportOut << endl << "Quadratic calibration residuals...." << endl;
      ReportOut << "Centroid (ch)\t\tList Energy (keV)\t\tCalibration Energy (keV)\t\tResidual (keV)" << endl;
      NumFits = 0;
      for (i = 0; i <= Config.Sources[Settings.Source].size(); i++) {
         if (Fit->FitSuccess[i] == 1) {
            CalibEn = Fit->QuadGainFit[0] + (Fit->QuadGainFit[1] * (Fit->PeakFits[i].Mean / Settings.Integration)) +
                (pow((Fit->PeakFits[i].Mean / Settings.Integration), 2) * Fit->QuadGainFit[2]);
            ReportOut << Fit->PeakFits[i].Mean << "\t\t\t" << Fit->PeakFits[i].Energy << "\t\t\t";
            ReportOut << CalibEn << "\t\t\t" << CalibEn - Fit->PeakFits[i].Energy << endl;
            Residuals[NumFits] = CalibEn - Fit->PeakFits[i].Energy;
            NumFits++;
         }
      }
      // Residual from linear fit
      NumFits = 0;
      ReportOut << endl << "Linear calibration residuals...." << endl;
      ReportOut << "Centroid (ch)\t\tList Energy (keV)\t\tCalibration Energy (keV)\t\tResidual (keV)" << endl;
      for (i = 0; i <= Config.Sources[Settings.Source].size(); i++) {
         if (Fit->FitSuccess[i] == 1) {
            CalibEn = Fit->LinGainFit[0] + (Fit->LinGainFit[1] * (Fit->PeakFits[i].Mean / Settings.Integration));
            ReportOut << Fit->PeakFits[i].Mean << "\t\t\t" << Fit->PeakFits[i].Energy << "\t\t\t";
            ReportOut << CalibEn << "\t\t\t" << CalibEn - Fit->PeakFits[i].Energy << endl;
            Energies[NumFits] = Fit->PeakFits[i].Energy;
         }
      }
      TGraphErrors CalibResidual(NumFits, Energies, Residuals);
      if (PLOT_RESIDUAL) {
         cCalib1a->cd(2);
         CalibResidual.SetMarkerColor(2);
         CalibResidual.SetMarkerStyle(20);
         CalibResidual.SetMarkerSize(1.0);
         CalibResidual.SetTitle("Residual from quadratic calibration");
         CalibResidual.Draw("AP");
         //CalibResidual.Draw();
         cCalib1a->Modified();
         cCalib1a->Update();
         App->Run(1);
         //cCalib1->cd(1);
      }
   }
   if(Config.CalFile) {
      
   }

}

int WriteCalFile(SpectrumFit *Fit, ofstream &CalFileOut, std::string HistName, FitSettings Settings) {



   return 0;
}
