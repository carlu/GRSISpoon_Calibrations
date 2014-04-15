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
//#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TFolder.h>


// TriScope libraries
//#include "TTigFragment.h"

// My libraries
#include "Calib.h"
#include "Main.h"
#include "CalibTools.h"

//Functions
// Convert crystal Name/number


// Fitting stuff
std::string FitOptions = ("RQE");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
std::string Opts;
//float InitialSigma = 0.0;
float InitialSigma = 50.0;
float InitialGain = 0.16;


extern TApplication *App;       // Pointer to root environment for plotting etc

extern TCanvas *cCalib1, *cCalib1a, *cCalib2, *cWave1, *ctemp;

int CalibSpectra(std::string filename)
{

   // Variables
   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   int Source = 0;
   int FitSuccess = 0;
   bool PlotOn = 0;
   int FileType = 0;
   int Integration = 0;
   float Dispersion = 0.0;
   int i, x;
   char CharBuf[CHAR_BUFFER_SIZE];
   string HistName;
   string OutputName;
   char Charge[10];
   char Colours[] = "BGRW";
   TH1F *GainPlot, *OffsetPlot, *QuadPlot;
   TH1F *GainHist, *OffsetHist, *QuadHist;
   int ItemNum = 0;
   ofstream GainOut;
   ofstream ReportOut;
   ofstream WaveOut;
   ofstream WaveReportOut;
   ofstream CalFileOut;

   int NumFits;

   float CalibEn = 0.0;
   string TestString = "his";

   bool FileTypeFound = 0;

   // Check file type
   if (strncmp(filename.c_str(), Config.CalOut.c_str(), 8) == 0) {
      FileType = 1;             // File is from output of Calib()
      FileTypeFound = 1;
   }
   if (strncmp(filename.c_str(), TestString.c_str(), 3) == 0) {
      FileType = 2;             // file is from TIGRESS DAQ analyser
      Config.CalWave = 0;
      FileTypeFound = 1;
   }

   if (FileTypeFound == 0) {
      cout << "Histogram file type not determined!" << endl;
      return -1;
   }
   // File
   TFile *file = TFile::Open(filename.c_str());
   if (file->IsOpen()) {
      if (Config.PrintBasic) {
         cout << filename << " opened!" << endl;
      }
   } else {
      if (Config.PrintBasic) {
         cout << "Failed to open " << filename << "!" << endl;
      }
      return -1;
   }

   // Open TFolder if required
   TFolder *Folder;
   if (FileType == 2) {
      Folder = (TFolder *) file->FindObjectAny("histos");
   }
   if (!Folder && FileType == 2) {
      cout << "Folder not loaded" << endl;
      return -1;
   }
   // Prepare output root file.
   std::string tempstring;
   TFile *outfile;
   if (Config.WriteFits) {
      tempstring = Config.OutPath + Config.CalSpecOut;
      outfile = TFile::Open(tempstring.c_str(), "RECREATE");
   }
   // Prepare energy calibration output files
   if (Config.CalEnergy) {
      tempstring = Config.OutPath + "GainsOut.txt";
      GainOut.open(tempstring.c_str());
      if (Config.CalReport) {
         tempstring = Config.OutPath + "CalibrationReport.txt";
         ReportOut.open(tempstring.c_str());
      }
   }
   // Prepart Wave calibration output files
   if (Config.CalWave) {
      tempstring = Config.OutPath + "WaveGainsOut.txt";
      WaveOut.open(tempstring.c_str());
      if (Config.CalReport) {
         tempstring = Config.OutPath + "WaveCalibrationReport.txt";
         WaveReportOut.open(tempstring.c_str());
      }
   }
   // Prepare .cal file
   if (Config.CalFile) {
      tempstring = Config.OutPath + "EnergyCalibration.cal";
      CalFileOut.open(tempstring.c_str());
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

   GainHist = new TH1F("Gain Histogram", "Histogram of Gain of all fitted channels", 512, 0, 0.5);
   GainHist->GetXaxis()->SetTitle("keV/ch");
   OffsetHist = new TH1F("Offset Histogram", "Histogram of Offset of all fitted channels", 512, -5, 5);
   OffsetHist->GetXaxis()->SetTitle("keV");
   QuadHist =
       new TH1F("Quadratic Histogram", "Histogram of Quadratic component of all fitted channels", 512, -0.000001,
                0.000001);
   QuadHist->GetXaxis()->SetTitle("keV/ch^2");
   TH1F *Histo;
   // If we're calibrating the FPGA energy...
   if (Config.CalEnergy) {
      // Loop all gamma detectors 
      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            for (Seg = 0; Seg <= SEGS + 1; Seg++) {

               if (Config.CalList[Clover - 1][Crystal][Seg] == 0) {     // check if this channel is to be fitted.
                  continue;
               }
               // Set source and histogram name and output name based on seg number and file type
               if (FileType == 1) {
                  switch (Seg) {
                  case 0:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da Chg", Clover, Num2Col(Crystal), Seg);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     Source = Config.SourceNumCore;
                     break;
                  case 9:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db Chg", Clover, Num2Col(Crystal), 0);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     Source = Config.SourceNumCore;
                     break;
                  default:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx Chg", Clover, Num2Col(Crystal), Seg);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     if (Seg < 5) {
                        Source = Config.SourceNumFront;
                     } else {
                        Source = Config.SourceNumBack;
                     }
                     break;
                  }
               }
               if (FileType == 2) {
                  switch (Crystal) {
                  case 0:
                     x = 0;
                     break;
                  case 1:
                     x = 20;
                     break;
                  case 2:
                     x = 30;
                     break;
                  case 3:
                     x = 50;
                     break;
                  }
                  switch (Seg) {
                  case 0:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "Chrg0%03d", ((Clover - 1) * 60) + x);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     OutputName += "\t";
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da Chg", Clover, Num2Col(Crystal), Seg);
                     OutputName += CharBuf;
                     Source = Config.SourceNumCore;
                     break;
                  case 9:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "Chrg0%03d", ((Clover - 1) * 60) + x + 9);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     OutputName += "\t";
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db Chg", Clover, Num2Col(Crystal), 0);
                     OutputName += CharBuf;
                     Source = Config.SourceNumCore;
                     break;
                  default:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "Chrg0%03d", ((Clover - 1) * 60) + x + Seg);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     OutputName += "\t";
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx Chg", Clover, Num2Col(Crystal), Seg);
                     OutputName += CharBuf;
                     if (Seg < 5) {
                        Source = Config.SourceNumFront;
                     } else {
                        Source = Config.SourceNumBack;
                     }
                     break;
                  }


               }
               // Load histogram
               if (FileType == 1) {     // Files where the histogram can be found directly
                  Histo = (TH1F *) file->FindObjectAny(HistName.c_str());
               } else {         // Files where it needs to be puled from a TFolder i.e. TIGRESS DAQ files
                  Histo = (TH1F *) Folder->FindObjectAny(HistName.c_str());
               }

               if (Histo) {
                  if (Config.PrintVerbose) {
                     cout << endl << "------------------------------------" << endl;
                     cout << "Hist " << HistName << " loaded" << endl;
                     cout << "------------------------------------" << endl << endl;
                  }
                  // Check if plot should be active for this channel
                  PlotOn = 0;
                  if (Config.PlotFits) {
                     if (Config.CalibPlots[Clover - 1][Crystal][Seg]) {
                        PlotOn = 1;
                     }
                  }
                  // Perform Fit                  
                  SpectrumFit Fit = { 0 };
                  FitSettings Settings = { 0 };

                  Settings.Source = Source;
                  if (FileType == 1) {
                     Settings.Integration = INTEGRATION;
                     Settings.Dispersion = float (CHARGE_BINS) / float (CHARGE_MAX);
                  }
                  if (FileType == 2) {  // histogram file from analyser
                     Settings.Integration = 1;
                     Settings.Dispersion = 1;
                  }
                  Settings.SearchSigma = EN_SEARCH_SIGMA;
                  Settings.SearchThresh = EN_SEARCH_THRESH;
                  Settings.SigmaEstZero = ENERGY_SIGMA_ZERO;
                  Settings.SigmaEst1MeV = ENERGY_SIGMA_1MEV;
                  Settings.FitZero = Config.FitZero;
                  Settings.PlotOn = PlotOn;

                  FitSuccess = FitGammaSpectrum(Histo, &Fit, Settings);

                  // If fit succesful, generate output....
                  if (FitSuccess > 0) {
                     switch (Crystal) { // Calculate channel number (old TIGRESS DAQ numbering)
                     case 0:
                        ItemNum = ((Clover - 1) * 60) + Seg;
                        break;
                     case 1:
                        ItemNum = ((Clover - 1) * 60) + 20 + Seg;
                        break;
                     case 2:
                        ItemNum = ((Clover - 1) * 60) + 30 + Seg;
                        break;
                     case 3:
                        ItemNum = ((Clover - 1) * 60) + 50 + Seg;
                        break;
                     default:
                        ItemNum = 1000;
                        break;
                     }
                     GainPlot->SetBinContent(ItemNum + 1, Fit.QuadGainFit[1]);  // ItemNum+1 to skip bin 0 which is underflow
                     OffsetPlot->SetBinContent(ItemNum + 1, Fit.QuadGainFit[0]);
                     QuadPlot->SetBinContent(ItemNum + 1, Fit.QuadGainFit[2]);

                     GainHist->Fill(Fit.QuadGainFit[1]);
                     OffsetHist->Fill(Fit.QuadGainFit[0]);
                     QuadHist->Fill(Fit.QuadGainFit[2]);
                  }
                  // Now print reports on results of fits and calibration.
                  if (FitSuccess > 0) {
                     if (FitSuccess < 3 || FORCE_LINEAR) {
                        GainOut << OutputName << ":\t" << Fit.LinGainFit[0];
                        GainOut << "\t" << Fit.LinGainFit[1] << endl;
                     } else {
                        GainOut << OutputName << ":\t" << Fit.QuadGainFit[0] << "\t";
                        GainOut << Fit.QuadGainFit[1] << "\t" << Fit.QuadGainFit[2] << endl;
                     }
                  } else {
                     //GainOut << HistName << " Fail!!!" << endl;
                  }
                  // Write full calibration report
                  if (FitSuccess > 0) {
                     CalibrationReport(&Fit, ReportOut, OutputName, Settings);
                  } else {
                     ReportOut << endl << "------------------------------------------" << endl << OutputName << endl <<
                         "------------------------------------------" << endl << endl;
                     ReportOut << "Fail Fail Fail! The calibration has failed!" << endl;
                  }
                  // Write .cal file for GRSISpoon
                  if (Config.CalFile) {
                     WriteCalFile(&Fit, CalFileOut, HistName, Settings);
                  }
                  // Write histo with fits
                  if (Config.WriteFits) {
                     outfile->cd();
                     Histo->Write();
                  }
               } else {
                  if (Config.PrintBasic) {
                     cout << endl << "Hist " << HistName << " failed to load." << endl;
                  }
               }
            }
         }
      }
   }
   // Now run the fit for the waveform spectrum if required
   if (Config.CalWave) {
      if (Config.PrintVerbose) {
         if (Config.PrintBasic) {
            cout << "-------------------" << endl << "Now fitting Wave Energy Spectra" << endl << "-------------------"
                << endl;
         }
      }
      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            for (Seg = 0; Seg <= SEGS + 1; Seg++) {

               if (Config.CalList[Clover - 1][Crystal][Seg] == 0) {     // check if this channel is to be fitted.
                  continue;
               }

               if (Config.PrintVerbose) {
                  cout << endl << "--------------------------------------" << endl;
                  cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                  cout << "--------------------------------------" << endl;
               }
               if (Seg == 0 || Seg == 9) {
                  Source = SOURCE_NUM_CORE;
                  if (Seg == 0) {
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da WaveChg", Clover, Colours[Crystal], Seg);
                     HistName = CharBuf;
                  } else {
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db WaveChg", Clover, Colours[Crystal], 0);
                     HistName = CharBuf;
                  }
               } else {
                  snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx WaveChg", Clover, Colours[Crystal], Seg);
                  HistName = CharBuf;
                  if (Seg < 5) {
                     Source = SOURCE_NUM_FRONT;
                  } else {
                     Source = SOURCE_NUM_BACK;
                  }
               }

               TH1F *Histo = (TH1F *) file->FindObjectAny(HistName.c_str());
               // Check if plot should be active for this channel
               PlotOn = 0;
               if (Config.PlotFits) {
                  if (Config.CalibPlots[Clover - 1][Crystal][Seg]) {
                     PlotOn = 1;
                  }
               }
               // Perform Fit                  
               SpectrumFit WaveFit = { 0 };
               FitSettings Settings = { 0 };

               Settings.Source = Source;
               Settings.Integration = 1;
               Settings.Dispersion = float (CHARGE_BINS) / float (WAVE_CHARGE_MAX);
               Settings.SearchSigma = WAVE_SEARCH_SIGMA;
               Settings.SearchThresh = WAVE_SEARCH_THRESH;
               Settings.SigmaEstZero = WAVE_SIGMA_ZERO;
               Settings.SigmaEst1MeV = WAVE_SIGMA_1MEV;
               Settings.FitZero = Config.FitZero;
               Settings.PlotOn = PlotOn;

               FitSuccess = FitGammaSpectrum(Histo, &WaveFit, Settings);

               if (FitSuccess > 0) {
                  if (FitSuccess < 3 || FORCE_LINEAR) {
                     WaveOut << HistName << "\t" << WaveFit.LinGainFit[0];
                     WaveOut << "\t" << WaveFit.LinGainFit[1] << endl;
                  } else {
                     WaveOut << HistName << ":\t" << WaveFit.QuadGainFit[0] << "\t";
                     WaveOut << WaveFit.QuadGainFit[1] << "\t" << WaveFit.QuadGainFit[2] << endl;
                  }
               } else {
                  //WaveOut << HistName << " Fail!!!" << endl;
               }
               if (Config.CalReport) {
                  if (FitSuccess > 0) {
                     CalibrationReport(&WaveFit, WaveReportOut, HistName, Settings);
                  } else {
                     WaveReportOut << endl << "------------------------------------------" << endl << HistName << endl
                         << "------------------------------------------" << endl << endl;
                     WaveReportOut << "Fail Fail Fail! The calibration has failed!" << endl;
                  }
               }
               if (Config.WriteFits) {
                  outfile->cd();
                  Histo->Write();
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

   if (Config.WriteFits) {
      outfile->Close();
   }

   if (Config.CalEnergy) {
      GainOut.close();
   }
   if (Config.CalReport) {
      ReportOut.close();
   }

   return 0;

}
