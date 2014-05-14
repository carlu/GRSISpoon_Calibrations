// To  compile: g++ SortHistos.C CalibTools.C Options.C Utils.C -I$GRSISYS/include --std=c++0x -o SortHistos $GRSISYS/libraries/TigFormat/libFormat.so $GRSISYS/libraries/libCalManager.so $GRSISYS/libraries/libRootIOManager.so -O0 `root-config --cflags --libs`  -lSpectrum -lgsl -lgslcblas -g

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
#include "Options.h"
#include "SortTrees.h"
#include "Calib.h"
#include "CalibTools.h"
#include "Utils.h"

//Functions
int CalibrateFiles(); 
int FitHistoFile(TFile *file, int FileType, MasterFitMap FitMap, MasterFitMap WaveFitMap);
int CalibrateChannel(ChannelFitMap Fits);

// Fitting stuff
std::string FitOptions = ("RQE");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
std::string Opts;
//float InitialSigma = 0.0;
float InitialSigma = 50.0;
float InitialGain = 0.16;

TApplication *App;       // Pointer to root environment for plotting etc

// Canvas for plots
TCanvas *cCalib1, *cCalib1a, *cCalib2, *cWave1, *ctemp;

// Spectra for Calibration
TH1F *GainPlot, *OffsetPlot, *QuadPlot;
TH1F *GainHist, *OffsetHist, *QuadHist;

int main(int argc, char **argv)
{
   // Set default and read custom options
   LoadDefaultSettings();
   if (ReadCommandLineSettings(argc, argv) < 0) {
      cout << "Failed to configure the run - exiting!" << endl;
      return -1;
   }
   
   // Set options for histo stats
   gStyle->SetOptStat("iouRMen");

   // create root environment for interacting with plots etc
   App = new TApplication("Output", 0, NULL);

   // Timing
   TStopwatch StopWatch;
   StopWatch.Start();

   // Initialise TCanvas's
   if (Config.RunSpecCal == 1) {
      cCalib1 = new TCanvas("cCalib1", "Fit", 800, 600);        // Canvas for spectrum plots

      cCalib1a = new TCanvas("cCalib1a", "Calibration", 800, 600);      // Canvas for spectrum plots
      cCalib1a->Divide(1, 2);

      cCalib2 = new TCanvas("cCalib2", "Calibration Summary", 800, 600);        // Canvas for gain plots and histograms
      cCalib2->Divide(2, 3);
      cCalib2->Update();

      cCalib1->cd();
   }
   
   // Check what we are supposed to be doing call function
   if (Config.RunSpecCal == 1) {
      if(CalibrateFiles() > 0) {
         cout << "CalibrateFiles() failed." << endl;
         return 1;
      }
      //filename = Config.files.at(i);
   }
   
   if(Config.RunSpecEff == 1)  {
      
   }  
   return 0;
}



int CalibrateFiles() {
   
   // Variables, Constants, etc
   int i, j, x;
   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   unsigned int ChainEvent, TreeEvent;
   
   int FileType;
   int CalibSuccess;
   
   char Charge[10];
   char Colours[] = "BGRW";
   int ItemNum = 0;
   ofstream GainOut;
   ofstream ReportOut;
   ofstream WaveOut;
   ofstream WaveReportOut;
   ofstream CalFileOut;
   int NumFits;
   float CalibEn = 0.0;
   bool FileTypeFound = 0;
   
   string filename, runname;
   
   std::string tempstring;
   TFile *file;

   std::vector <int> ChanVector;  // vector to store (Clover,Crystal,Seg) for use as map key
   ChanVector.push_back(0); ChanVector.push_back(0); ChanVector.push_back(0);  // populate chanvector with zeroes.
   
   // Maps to store all the fits in.
   MasterFitMap FitMap;
   MasterFitMap WaveFitMap;
   
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
   
   // Check we have sources for all runs
   if(Config.files.size() != Config.SourceNumCore.size()) {
      cout << "For histogram calibration the number of files (" << Config.files.size() << ") must equal number of sources (" << Config.SourceNumCore.size()  << ")" << endl;
      return -1;
   }
   
   // loop input files
   for (i = 0; i < Config.files.size(); i++) {
   
      if (Config.PrintBasic) {
         cout << "Attempting offline calibration on histograms in file: " << Config.files.at(i) << endl;
      }
      
      filename = Config.files.at(i);
   
      // Check file type
      if (strncmp(filename.c_str(), Config.CalName.c_str(), Config.CalName.size()) == 0) {
         FileType = 1;             // File is from output of Calib()
         FileTypeFound = 1;
      }
      
      if (strncmp(filename.c_str(), Config.AnaName.c_str(), Config.AnaName.size()) == 0) {
         FileType = 2;             // file is from TIGRESS DAQ analyser
         Config.CalWave = 0;
         FileTypeFound = 1;
      }
      
      if (FileTypeFound == 0) {
         if (Config.PrintBasic) {
            cout << filename << ": Histogram file type not determined, skipping." << endl;
            cout << "(Note, using ./ at the begining of the filename screws this up.  Must fix.)" << endl;
         }
         continue;
      }
   
      // Input File
      file = TFile::Open(filename.c_str(),"READ");
      if (file->IsOpen()) {
         if (Config.PrintBasic) {
            cout << filename << " opened!" << endl;
         }
      } else {
         if (Config.PrintBasic) {
            cout << "Failed to open " << filename << "!" << endl;
         }
         continue;
      }
      
      // Now fit the histos in file
      FitHistoFile(file, FileType, FitMap, WaveFitMap);
   
   }
   
   // Now iterate over FitMap/WaveFitMap and perform calibrations
   // This should mayeb be a separte function to loop channels and another to actually calibrate.
   // Check fit was succesful before calling
   //CalibSuccess = CalibrateGammaSpectrum(&Fit, Settings);
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg <= SEGS + 1; Seg++) {
         
            if (Config.CalList[Clover - 1][Crystal][Seg] == 0) {     // check if this channel is to be fitted.
               continue;
            }
            
            // Set Cl,Cr,Seg vector for use as key in Master fit map
            ChanVector.at(0) = Clover;
            ChanVector.at(1) = Crystal;
            ChanVector.at(2) = Seg;
            
            CalibSuccess = CalibrateChannel(FitMap[ChanVector]);
         
         }
      }   
   }
   
}

int FitHistoFile(TFile *file, int FileType, MasterFitMap FitMap, MasterFitMap WaveFitMap) {

   int i, j, x;
   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   int Source = 0;
   int Line;
   
   int SourceNumCore, SourceNumFront, SourceNumBack;
   int FitSuccess = 0;
   int CalibSuccess = 0;
   bool PlotOn = 0;
   bool PeakSelect=0;
   int Integration = 0;
   float Dispersion = 0.0;
   
   std::string tempstring;
   char CharBuf[CHAR_BUFFER_SIZE];
   string HistName;
   string OutputName;
   
   std::vector <int> ChanVector;  // vector to store (Clover,Crystal,Seg) for use as map key
   ChanVector.push_back(0); ChanVector.push_back(0); ChanVector.push_back(0);  // populate chanvector with zeroes.
   ChannelFitMap ChanFits;
   
   // ROOT objects
   TFile *outfile;
   TDirectory *dCharge, *dWaveCharge, *dCalibration = { 0 };
   TH1F *Histo;
   TFolder *Folder = 0;
   
   // Open TFolder if required
   if (FileType == 2) {
      Folder = (TFolder *) file->FindObjectAny("histos");
   }
   if (!Folder && FileType == 2) {
      cout << "Folder not loaded" << endl;
      return -1;
   }
   // Prepare output root file.
   if (Config.WriteFits) {
      tempstring = Config.OutPath + Config.CalSpecOut;
      outfile = TFile::Open(tempstring.c_str(), "RECREATE");
      dCharge = outfile->mkdir("Charge");
      dWaveCharge = outfile->mkdir("WaveCharge");
      dCalibration = outfile->mkdir("Calibration");
   }
   
   // Get source definition
   SourceNumCore = Config.SourceNumCore.at(i);
   SourceNumFront= Config.SourceNumFront.at(i);
   SourceNumBack = Config.SourceNumBack.at(i);
   
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
               // This should be factored out to a separate function
               if (FileType == 1) {
                  switch (Seg) {
                  case 0:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da Chg", Clover, Num2Col(Crystal), Seg);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     Source = SourceNumCore;
                     break;
                  case 9:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db Chg", Clover, Num2Col(Crystal), 0);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     Source = SourceNumCore;
                     break;
                  default:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx Chg", Clover, Num2Col(Crystal), Seg);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     if (Seg < 5) {
                        Source = SourceNumFront;
                     } else {
                        Source = SourceNumBack;
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
                     Source = SourceNumCore;
                     break;
                  case 9:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "Chrg0%03d", ((Clover - 1) * 60) + x + 9);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     OutputName += "\t";
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db Chg", Clover, Num2Col(Crystal), 0);
                     OutputName += CharBuf;
                     Source = SourceNumCore;
                     break;
                  default:
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "Chrg0%03d", ((Clover - 1) * 60) + x + Seg);
                     HistName = CharBuf;
                     OutputName = CharBuf;
                     OutputName += "\t";
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx Chg", Clover, Num2Col(Crystal), Seg);
                     OutputName += CharBuf;
                     if (Seg < 5) {
                        Source = SourceNumFront;
                     } else {
                        Source = SourceNumBack;
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
                  // Check if Manual peak select should be active
                  if (Config.ManualPeakSelect[Clover - 1][Crystal][Seg]) {
                     PeakSelect  = 1;
                  }
                  else {
                     PeakSelect = 0;
                  }
                  
                  // Initialise FitSettings and HistoFit               
                  HistoFit Fit;
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
                  Settings.PeakSelect = PeakSelect;
                  Settings.BackupPeakSelect = 0;
                  
                  // Perform fit
                  FitSuccess = FitGammaSpectrum(Histo, &Fit, Settings);
                  
                  // Write fit results to map
                  // First add all fitted lines to ChannelFit map
                  ChanFits.clear();
                  for(Line=0; Line<Fit.PeakFits.size(); Line ++) {
                     ChanFits.insert(ChannelFitPair(Fit.PeakFits.at(Line).Energy,Fit.PeakFits.at(Line)));
                  }
                  // Set Cl,Cr,Seg vector for use as key in Master fit map
                  ChanVector.at(0) = Clover;
                  ChanVector.at(1) = Crystal;
                  ChanVector.at(2) = Seg;
                  // Check if map entry exists for this channel and insert if not.
                  if(!FitMap.count(ChanVector)) {  // If entry for this channel does not exist...
                     FitMap.insert(MasterFitPair(ChanVector,ChanFits));  // create it....
                  }
                  else {                                                   // else....
                     
                     FitMap[ChanVector].insert(ChanFits.begin(),ChanFits.end());      // add to it.
                  }
                  
                  // Write histo with fits
                  if (Config.WriteFits) {
                     dCharge->cd();
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

               if (Seg == 0 || Seg == 9) {
                  Source = SourceNumCore;
                  if (Seg == 0) {
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da WaveChg", Clover, Num2Col(Crystal), Seg);
                     HistName = CharBuf;
                  } else {
                     snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db WaveChg", Clover, Num2Col(Crystal), 0);
                     HistName = CharBuf;
                  }
               } else {
                  snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx WaveChg", Clover, Num2Col(Crystal), Seg);
                  HistName = CharBuf;
                  if (Seg < 5) {
                     Source = SourceNumFront;
                  } else {
                     Source = SourceNumBack;
                  }
               }

               TH1F *Histo = (TH1F *) file->FindObjectAny(HistName.c_str());
               
               if(Histo) {
               
                  if (Config.PrintVerbose) {
                     cout << endl << "--------------------------------------" << endl;
                     cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                     cout << "--------------------------------------" << endl;
                  }
                  // Check if plot should be active for this channel
                  PlotOn = 0;
                  if (Config.PlotFits) {
                     if (Config.CalibPlots[Clover - 1][Crystal][Seg]) {
                        PlotOn = 1;
                     }
                  }
                  // Check if Manual peak select should be active
                  if (Config.ManualPeakSelect[Clover - 1][Crystal][Seg]) {
                     PeakSelect  = 1;
                  }
                  else {
                     PeakSelect = 0;
                  }
                     
                  // Perform Fit                  
                  HistoFit WaveFit;
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
                  Settings.PeakSelect = PeakSelect;
                  Settings.BackupPeakSelect = 0;
                  
                  FitSuccess = FitGammaSpectrum(Histo, &WaveFit, Settings);
                  
                  // Write fit results to map
                  // First add all fitted lines to ChannelFit map
                  ChanFits.clear();
                  for(Line=0; Line<WaveFit.PeakFits.size(); Line ++) {
                     ChanFits.insert(ChannelFitPair(WaveFit.PeakFits.at(Line).Energy,WaveFit.PeakFits.at(Line)));
                  }
                  // Set Cl,Cr,Seg vector for use as key in Master fit map
                  ChanVector.at(0) = Clover;
                  ChanVector.at(1) = Crystal;
                  ChanVector.at(2) = Seg;
                  // Check if map entry exists for this channel and insert if not.
                  if(!FitMap.count(ChanVector)) {  // If entry for this channel does not exist...
                     FitMap.insert(MasterFitPair(ChanVector,ChanFits));  // create it....
                  }
                  else {                                                   // else....
                     
                     FitMap[ChanVector].insert(ChanFits.begin(),ChanFits.end());      // add to it.
                  }
                  
                  if (Config.WriteFits) {
                     dWaveCharge->cd();
                     Histo->Write();
                  }
               }
            }
         }
      }

   }



   if (Config.WriteFits) {
      outfile->Close();
   }

   

   return 0;
   
}



int CalibrateChannel(ChannelFitMap Fits) {

   int Clover, Crystal, Seg;
   int ItemNum;

   // If fit succesful, generate output....
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
   
   // Now print reports on results of fits and calibration.
   if (Fit.LinesUsed > 1) {
      if (Fit.LinesUsed < 3 || FORCE_LINEAR) {
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
   if (Fit.LinesUsed > 0) {
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
   
   
   if (Config.CalEnergy) {
      GainOut.close();
   }
   if (Config.CalReport) {
      ReportOut.close();
   }

}

