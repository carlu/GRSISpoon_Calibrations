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
#include "CalibTools.h"
//#include "FitGammaSpectrumWave.h"

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


TApplication *App;              // Pointer to root environment for plotting etc

TCanvas *cCalib1, *cCalib2;     //, *cCalib3;

int main(int argc, char **argv)
{

   // Variables
   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   int Source = 0;
   int FitSuccess = 0;
   int PlotOn = 0;
   int Integration = 0;
   int i;
   char HistName[1024];
   char Charge[10];
   char Colours[] = "BGRW";
   TH1F *GainPlot, *OffsetPlot, *QuadPlot;
   TH1F *GainHist, *OffsetHist, *QuadHist;
   int ItemNum = 0;
   ofstream GainOut;
   ofstream ReportOut;
   ofstream WaveOut;

   int NumFits;
   float Energies[NUM_LINES], Residuals[NUM_LINES];

   float CalibEn = 0.0;

   gStyle->SetOptStat("iouRMen");

   // Set up plots
   App = new TApplication("Output", 0, NULL);   // creates root environment for interacting with plots etc
   if (PLOT_FITS || PLOT_CALIB || PLOT_CALIB_SUMMARY || PLOT_RESIDUAL) {
      cCalib1 = new TCanvas("cCalib1", "Fit and Calibration", 800, 600);        // Canvas for spectrum plots
      cCalib1->Divide(1, 3);

      cCalib2 = new TCanvas("cCalib2", "Calibration Summary", 800, 600);        // Canvas for gain plots and histograms
      cCalib2->Divide(2, 3);
      cCalib2->Update();

      //cCalib3 = new TCanvas("cCalib3", "Calibration Residual", 800, 600);  // Canvas for residual plit
      //cCalib3->Update();   

      cCalib1->cd(1);
   }
   // Files
   if (argc < 2) {
      cout << "Filename!" << endl;
      return 0;
   }
   
   TFile *file = TFile::Open(argv[1]);
   if(file->IsOpen()) {
      cout << argv[1] << " opened!" << endl;
   }
   else {
      cout << "Failed to open " << argv[1] << "!" << endl;
      return 0;
   }

   if (OUTPUT_GAIN) {
      GainOut.open("GainsOut.txt");
   }
   if (OUTPUT_REPORT) {
      ReportOut.open("CalibrationReport.txt");
   }
   if (FIT_WAVE_EN) {
      WaveOut.open("WaveGainsOut.txt");
   }
   // Histograms
   GainPlot = new TH1F("Gains", "Gain of all fitted channels", 1001, -0.5, 1000.5);
   GainPlot->GetYaxis()->SetTitle("keV/ch");
   GainPlot->GetXaxis()->SetTitle("Channel");
   OffsetPlot = new TH1F("Offsets", "Offset of all fitted channels", 1001, -0.5, 1000.5);
   OffsetPlot->GetYaxis()->SetTitle("keV");
   OffsetPlot->GetXaxis()->SetTitle("Channel");
   QuadPlot = new TH1F("Quads", "Quadratic component of all fitted channels", 1001, -0.5, 1000.5);
   QuadPlot->GetYaxis()->SetTitle("keV/ch^2");
   QuadPlot->GetXaxis()->SetTitle("Channel");

   GainHist = new TH1F("Gain Histogram", "Histogram of Gain of all fitted channels", 512, 0, 0.3);
   GainHist->GetXaxis()->SetTitle("keV/ch");
   OffsetHist = new TH1F("Offset Histogram", "Histogram of Offset of all fitted channels", 512, -5, 5);
   OffsetHist->GetXaxis()->SetTitle("keV");
   QuadHist =
       new TH1F("Quadratic Histogram", "Histogram of Quadratic component of all fitted channels", 512, -0.000001,
                0.000001);
   QuadHist->GetXaxis()->SetTitle("keV/ch^2");

   // If we're calibrating the FPGA energy...
   if (FIT_EN) {
      // Loop all gamma detectors 
      //for (Clover = 0; Clover < CLOVERS; Clover++) {
      for (Clover = 11; Clover < 12; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            for (Seg = 0; Seg <= SEGS + 1; Seg++) {

               // Set source and histogram name based on seg number
               switch (Seg) {
               case 0:
                  sprintf(HistName, "TIG%02d%cN%02da Chg", Clover + 1, Colours[Crystal], Seg);
                  Source = SOURCE_NUM_CORE;
                  break;
               case 9:
                  sprintf(HistName, "TIG%02d%cN%02db Chg", Clover + 1, Colours[Crystal], 0);
                  Source = SOURCE_NUM_CORE;
                  break;
               default:
                  sprintf(HistName, "TIG%02d%cP%02dx Chg", Clover + 1, Colours[Crystal], Seg);
                  if (Seg < 5) {
                     Source = SOURCE_NUM_FRONT;
                  } else {
                     Source = SOURCE_NUM_BACK;
                  }
                  break;
               }

               // Load histogram
               TH1F *Histo = (TH1F *) file->FindObjectAny(HistName);
               if (Histo) {
                  cout << endl << "------------------------------------" << endl;
                  cout << "Hist " << HistName << " loaded" << endl;

                  // Check if plot should be active for this channel
                  PlotOn = 0;
                  if (PLOT_FITS) {
                     if (PLOT_CLOVER == 0 || (Clover + 1) == PLOT_CLOVER) {
                        if (PLOT_CRYSTAL == 0 || (Crystal + 1) == PLOT_CRYSTAL) {
                           if (PLOT_SEG == 0 || Seg == PLOT_SEG) {
                              PlotOn = 1;
                           }
                        }
                     }
                  }
                  
                  // Perform Fit
                  SpectrumFit Fit = {0};
                  Integration = INTEGRATION;
                  FitSuccess = FitGammaSpectrum(Histo, &Fit, Source, Integration, PlotOn);

                  // If fit succesful, generate output....
                  if (FitSuccess > 0) {
                     switch (Crystal) {  // Calculate channel number (old TIGRESS DAQ numbering)
                     case 0:
                        ItemNum = (Clover * 60) + Seg;
                        break;
                     case 1:
                        ItemNum = (Clover * 60) + 20 + Seg;
                        break;
                     case 2:
                        ItemNum = (Clover * 60) + 30 + Seg;
                        break;
                     case 3:
                        ItemNum = (Clover * 60) + 50 + Seg;
                        break;
                     default:
                        ItemNum = 1000;
                        break;
                     }
                     //cout << "C,C,S = " << Clover << ", " << Crystal << ", " << Seg << "  ItemNum = " << ItemNum << endl;
                     GainPlot->SetBinContent(ItemNum + 1, Fit.QuadGainFit[1]);  // ItemNum+1 to skip bin 0 which is underflow
                     OffsetPlot->SetBinContent(ItemNum + 1, Fit.QuadGainFit[0]);
                     QuadPlot->SetBinContent(ItemNum + 1, Fit.QuadGainFit[2]);

                     GainHist->Fill(Fit.QuadGainFit[1]);
                     OffsetHist->Fill(Fit.QuadGainFit[0]);
                     QuadHist->Fill(Fit.QuadGainFit[2]);
                  }
                  // Now print reports on results of fits and calibration.
                  if (OUTPUT_GAIN) {
                     if (FitSuccess > 0) {
                        //GainOut << HistName << " " << Fit.LinGainFit[0] << " +/- " << Fit.dLinGainFit[0];
                        //GainOut << " " << Fit.LinGainFit[1] << " +/- " << Fit.dLinGainFit[1] << " " <<   Fit.LinGainFit[2] << endl;
                        GainOut << HistName << " " << Fit.QuadGainFit[0] << " +/- " << Fit.dQuadGainFit[0] << " ";
                        GainOut << Fit.QuadGainFit[1] << " +/- " << Fit.
                            dQuadGainFit[1] << " " << Fit.QuadGainFit[2] << " +/- " << Fit.dQuadGainFit[2] << endl;
                     } else {
                        GainOut << HistName << " Fail!!!" << endl;
                     }
                  }
                  if (OUTPUT_REPORT) {
                     if (FitSuccess > 0) {
                        // Heading for this channel
                        ReportOut << endl << "------------------------------------------" << endl << HistName << endl <<
                            "------------------------------------------" << endl << endl;
                        // Peak fit info
                        for (i = 0; i < NUM_LINES; i++) {
                           if (Fit.FitSuccess[i] == 1) {
                              ReportOut << Fit.PeakFits[i].Energy << " keV" << endl;
                              ReportOut << Fit.PeakFits[i].Const << " +/- " << Fit.PeakFits[i].dConst << "\t";
                              ReportOut << Fit.PeakFits[i].Mean << " +/- " << Fit.PeakFits[i].dMean << "\t";
                              ReportOut << Fit.PeakFits[i].Const << " +/- " << Fit.PeakFits[i].dConst << endl;
                              ReportOut << Fit.PeakFits[i].ChiSq << "\t" << Fit.PeakFits[i].
                                  NDF << "\t" << Fit.PeakFits[i].ChiSq / Fit.PeakFits[i].NDF;
                              ReportOut << endl << endl;
                           }
                        }
                        // Calibration...
                        ReportOut << "Linear Solution: Offset = " << Fit.
                            LinGain[0] << " +/- " << Fit.dLinGain[0] << "\t";
                        ReportOut << "Gain = " << Fit.LinGain[1] << " +/- " << Fit.dLinGain[1] << endl;
                        ReportOut << "Linear Fit: Offset = " << Fit.
                            LinGainFit[0] << " +/- " << Fit.dLinGainFit[0] << "\t";
                        ReportOut << "Gain = " << Fit.LinGainFit[1] << " +/- " << Fit.dLinGainFit[1] << "\t";
                        ReportOut << "CSPD = " << Fit.LinGainFit[2] << endl;
                        ReportOut << "Quadratic Fit: Offset = " << Fit.
                            QuadGainFit[0] << " +/- " << Fit.dQuadGainFit[0] << "\t";
                        ReportOut << "Gain = " << Fit.QuadGainFit[1] << " +/- " << Fit.dQuadGainFit[1] << "\t";
                        ReportOut << "Quad = " << Fit.QuadGainFit[2] << " +/- " << Fit.dQuadGainFit[2] << "\t";
                        ReportOut << "CSPD = " << Fit.LinGainFit[3] << endl;

                        // Residual from quadratic fit
                        ReportOut << endl << "Check quadratic calibration...." << endl;
                        ReportOut << "Centroid (ch)\t\tList Energy (keV)\t\tCalibration Energy (keV)\t\tResidual (keV)"
                            << endl;
                        for (i = 0; i < NUM_LINES; i++) {
                           if (Fit.FitSuccess[i] == 1) {
                              CalibEn =
                                  Fit.QuadGainFit[0] + (Fit.QuadGainFit[1] * (Fit.PeakFits[i].Mean / INTEGRATION)) +
                                  (pow((Fit.PeakFits[i].Mean / INTEGRATION), 2) * Fit.QuadGainFit[2]);
                              ReportOut << Fit.PeakFits[i].Mean << "\t\t\t" << Fit.PeakFits[i].Energy << "\t\t\t";
                              ReportOut << CalibEn << "\t\t\t" << CalibEn - Fit.PeakFits[i].Energy << endl;
                           }
                        }

                        // Residual from linear fit
                        NumFits = 0;
                        ReportOut << endl << "Check linear calibration...." << endl;
                        ReportOut << "Centroid (ch)\t\tList Energy (keV)\t\tCalibration Energy (keV)\t\tResidual (keV)"
                            << endl;
                        for (i = 0; i < NUM_LINES; i++) {
                           if (Fit.FitSuccess[i] == 1) {
                              CalibEn = Fit.LinGainFit[0] + (Fit.LinGainFit[1] * (Fit.PeakFits[i].Mean / INTEGRATION));
                              ReportOut << Fit.PeakFits[i].Mean << "\t\t\t" << Fit.PeakFits[i].Energy << "\t\t\t";
                              ReportOut << CalibEn << "\t\t\t" << CalibEn - Fit.PeakFits[i].Energy << endl;
                              Energies[NumFits] = Fit.PeakFits[i].Energy;
                              Residuals[NumFits] = CalibEn - Fit.PeakFits[i].Energy;
                              NumFits++;
                           }
                        }
                        TGraphErrors CalibResidual(NumFits, Energies, Residuals);
                        if (PLOT_RESIDUAL) {
                           cCalib1->cd(3);
                           CalibResidual.SetMarkerColor(2);
                           CalibResidual.SetMarkerStyle(20);
                           CalibResidual.SetMarkerSize(1.0);
                           CalibResidual.SetTitle("Residual from linear calibration");
                           CalibResidual.Draw("AP");
                           //CalibResidual.Draw();
                           cCalib1->Modified();
                           cCalib1->Update();
                           App->Run(1);
                           //cCalib1->cd(1);
                        }
                     } else {
                        ReportOut << endl << "------------------------------------------" << endl << HistName << endl <<
                            "------------------------------------------" << endl << endl;
                        ReportOut << "Fail Fail Fail! The calibration has failed!" << endl;
                     }
                  }
               } else {
                  cout << endl << "Hist " << HistName << " failed to load." << endl;
               }
            }
         }
      }

      // Now run the fit for the waveform spectrum if required
      if (FIT_WAVE_EN) {
         if (VERBOSE) {
            cout << "-------------------" << endl << "Now fitting Wave Energy Spectra" << endl << "-------------------"
                << endl;
         }
         //for (Clover = 0; Clover < CLOVERS; Clover++) {
         for (Clover = 11; Clover < 12; Clover++) {
            for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
               for (Seg = 0; Seg <= SEGS + 1; Seg++) {

                  PlotOn = 1;
                  if (Seg == 0) {
                     PlotOn = 1;
                  }

                  if (VERBOSE) {
                     cout << endl << "--------------------------------------" << endl;
                     cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                     cout << "--------------------------------------" << endl;
                  }
                  if (Seg == 0 || Seg == 9) {
                     Source = SOURCE_NUM_CORE;
                     if (Seg == 0) {
                        sprintf(HistName, "TIG%02d%cN%02da WaveChg", Clover + 1, Colours[Crystal], Seg);
                     } else {
                        sprintf(HistName, "TIG%02d%cN%02db WaveChg", Clover + 1, Colours[Crystal], 0);
                     }
                  } else {
                     sprintf(HistName, "TIG%02d%cP%02dx WaveChg", Clover + 1, Colours[Crystal], Seg);
                     if (Seg < 5) {
                        Source = SOURCE_NUM_FRONT;
                     } else {
                        Source = SOURCE_NUM_BACK;
                     }
                  }

                  TH1F *Histo = (TH1F *) file->FindObjectAny(HistName);
                  // Perform Fit
                  SpectrumFit WaveFit = {};
                  Integration = 1;
                  FitSuccess = FitGammaSpectrum(Histo, &WaveFit, Source, Integration, PlotOn);

                  if (OUTPUT_GAIN) {
                     if (FitSuccess > 0) {
                        WaveOut << HistName << " " << WaveFit.LinGainFit[0] << " +/- " << WaveFit.dLinGainFit[0];
                        WaveOut << " " << WaveFit.LinGainFit[1] << " +/- " << WaveFit.dLinGainFit[1] << " " << WaveFit.LinGainFit[2] << endl;
                        WaveOut << HistName << " " << WaveFit.QuadGainFit[0] << " +/- " << WaveFit.dQuadGainFit[0] << " ";
                        WaveOut << WaveFit.QuadGainFit[1] << " +/- " << WaveFit.
                            dQuadGainFit[1] << " " << WaveFit.QuadGainFit[2] << " +/- " << WaveFit.dQuadGainFit[2] << endl;
                     } else {
                        WaveOut << HistName << " Fail!!!" << endl;
                     }
                  }
               }
            }
         }
      }
   } 

   if (PLOT_CALIB_SUMMARY) {
      //cCalib2->cd();
      cCalib2->cd(1);
      OffsetPlot->Draw();
      cCalib2->Modified();
      cCalib2->Update();
      //App->Run();
      cCalib2->cd(2);
      OffsetHist->Draw();
      cCalib2->Modified();
      cCalib2->Update();
      cCalib2->cd(3);
      GainPlot->Draw();
      cCalib2->Modified();
      cCalib2->Update();
      cCalib2->cd(4);
      GainHist->Draw();
      cCalib2->Modified();
      cCalib2->Update();
      cCalib2->cd(5);
      QuadPlot->Draw();
      cCalib2->Modified();
      cCalib2->Update();
      cCalib2->cd(6);
      QuadHist->Draw();
      cCalib2->Modified();
      cCalib2->Update();
      App->Run();
   }


   if (OUTPUT_GAIN) {
      GainOut.close();
   }
   if (OUTPUT_REPORT) {
      ReportOut.close();
   }

   return 0;

}
