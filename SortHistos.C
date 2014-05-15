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
int FitHistoFile(TFile *file, int FileType, int FileNum, MasterFitMap FitMap, MasterFitMap WaveFitMap);
int CalibrateChannel(ChannelFitMap Fits, FitSettings Settings);
int ConfigureEnergyFit(int Clover, int Crystal, int Seg, int FileType, int FileNum, FitSettings *Settings) ;
int ConfigureWaveEnFit(int Clover, int Crystal, int Seg,  int FileType, int FileNum, FitSettings *Settings) ;

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

// output streams 
ofstream GainOut;
ofstream ReportOut;
ofstream WaveOut;
ofstream WaveReportOut;
ofstream CalFileOut;

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
   int FileNum;
   int CalibSuccess;
   
   char Charge[10];
   char Colours[] = "BGRW";
   int ItemNum = 0;
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
   for (FileNum = 0; FileNum < Config.files.size(); FileNum++) {
   
      if (Config.PrintBasic) {
         cout << "Attempting offline calibration on histograms in file: " << Config.files.at(FileNum) << endl;
      }
      
      filename = Config.files.at(FileNum);
   
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
      FitHistoFile(file, FileType, FileNum, FitMap, WaveFitMap);
   
   }
   
   // Now iterate over FitMap/WaveFitMap and perform calibrations
   // This should mayeb be a separte function to loop channels and another to actually calibrate.
   // Check fit was succesful before calling
   //CalibSuccess = CalibrateGammaSpectrum(&Fit, Settings);
   
   if (Config.CalEnergy) {
      FitSettings Settings;
      FileNum = 0; // Reset as this value not needed for calibration, only fit
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
               
               ConfigureEnergyFit(Clover,Crystal, Seg, FileType, FileNum, &Settings);
               
               CalibSuccess = CalibrateChannel(FitMap[ChanVector], Settings);
            
            }
         }   
      }
   }
   
   if(Config.CalWave) {
      FitSettings Settings;
      FileNum = 0; // Reset as this value not needed for calibration, only fit
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
               
               ConfigureWaveEnFit(Clover,Crystal, Seg, FileType, FileNum, &Settings);
               
               CalibSuccess = CalibrateChannel(WaveFitMap[ChanVector], Settings);
            
            }
         }   
      }
   }
}

int FitHistoFile(TFile *file, int FileType, int FileNum, MasterFitMap FitMap, MasterFitMap WaveFitMap) {

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
   std::string HistName;
   std::string OutputName;
   
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
   
   // If we're calibrating the FPGA energy...
   if (Config.CalEnergy) {
      // Loop all gamma detectors 
      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            for (Seg = 0; Seg <= SEGS + 1; Seg++) {

               if (Config.CalList[Clover - 1][Crystal][Seg] == 0) {     // check if this channel is to be fitted.
                  continue;
               }
               
               FitSettings Settings = { 0 };
               
               ConfigureEnergyFit(Clover,Crystal, Seg, FileType, FileNum, &Settings);
               
               // Load histogram
               if (FileType == 1) {     // Files where the histogram can be found directly
                  Histo = (TH1F *) file->FindObjectAny(Settings.HistName.c_str());
               } else {         // Files where it needs to be puled from a TFolder i.e. TIGRESS DAQ files
                  Histo = (TH1F *) Folder->FindObjectAny(Settings.HistName.c_str());
               }

               if (Histo) {
                  if (Config.PrintVerbose) {
                     cout << endl << "------------------------------------" << endl;
                     cout << "Hist " << Settings.HistName << " loaded" << endl;
                     cout << "------------------------------------" << endl << endl;
                  }
                  
                  
                  
                  // Initialise FitSettings and HistoFit               
                  HistoFit Fit;
                  
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

               FitSettings Settings = { 0 };
               
               ConfigureWaveEnFit(Clover, Crystal, Seg, FileType, FileNum, &Settings);

               TH1F *Histo = (TH1F *) file->FindObjectAny(HistName.c_str());
               
               if(Histo) {
               
                  if (Config.PrintVerbose) {
                     cout << endl << "--------------------------------------" << endl;
                     cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                     cout << "--------------------------------------" << endl;
                  }
                  
                     
                  // Perform Fit                  
                  HistoFit WaveFit;
                  
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



int CalibrateChannel(ChannelFitMap Fits, FitSettings Settings) {
   
   // This function should copy the operation of the old CalibrateGammaSpectrum() in CalibTools.C
   // However it should be altered to use maps and vectors rather than arrays.
   int Clover, Crystal, Seg, i;
   int LinesUsed;
   int ItemNum;
   bool TestsPassed;
   
   float Energy = 0.0;
   float Energies[25], dEnergies[25];
   float Centroids[25], dCentroids[25];
   
   FitResult LineFit;
   HistoFit HistoCal;
   
   LinesUsed = 0;
   for (ChannelFitMapIt Line = Fits.begin(); Line != Fits.end(); Line++) {
      TestsPassed = 0;
      Energy = Line->first;
      LineFit = Line->second;
      
      HistoCal.PeakFits.push_back(LineFit); // Store fit for use in calibration and report 
      
      if (LineFit.Const > GAUS_HEIGHT_MIN) {  // Enough counts?
         if ((LineFit.ChiSq / LineFit.NDF) < GAUS_CSPD_MAX) {  // Good enouhg fit?
            // If passed, add this fit to calibration points
            Energies[LinesUsed] = Energy;
            dEnergies[LinesUsed] = 0.0; // ignore error here for now.
            Centroids[LinesUsed] = LineFit.Mean / Settings.Integration;
            dCentroids[LinesUsed] = LineFit.dMean / Settings.Integration;
            if (Config.PrintVerbose) {
               cout << Energies[LinesUsed] << " keV\t" << Centroids[LinesUsed] << " +/- " <<
                   dCentroids[LinesUsed] << " ch" << endl;
            }
            LinesUsed += 1;
            TestsPassed = 1;
            HistoCal.FitSuccess.push_back(1);
            if(LinesUsed > 23) {
               break;
            }
            
         }
      }
      if(TestsPassed==0) {
         HistoCal.FitSuccess.push_back(0);
      }
   }
   
   // If required, add point at 0ch=0keV
   if (Settings.FitZero) {
      Energies[LinesUsed] = 0.0;
      dEnergies[LinesUsed] = 0.0;
      Centroids[LinesUsed] = 0.0;
      dCentroids[LinesUsed] = ZERO_ERR;
      if (Config.PrintVerbose) {
         cout << Energies[LinesUsed] << " keV\t" << Centroids[LinesUsed] << " +/- " << dCentroids[LinesUsed] <<
             " ch" << endl;
      }
      LinesUsed += 1;
      FitResult ZeroFit;
      ZeroFit.Energy = 0.0;
      ZeroFit.Const = 0.0;
      ZeroFit.dConst = 0.0;
      ZeroFit.Mean = 0.0;
      ZeroFit.dMean = 1.0;
      ZeroFit.Sigma = 0.0;
      ZeroFit.dSigma = 0.0;
      ZeroFit.ChiSq = 0.0;
      ZeroFit. NDF = 1;
      HistoCal.PeakFits.push_back(ZeroFit);
      HistoCal.FitSuccess.push_back(1);
   }   
   
   HistoCal.LinesUsed = LinesUsed;
   
   // Now add points to a TGraph 
   TGraphErrors CalibPlot(LinesUsed, Centroids, Energies, dCentroids, dEnergies);
   if (Config.PrintVerbose) {
      cout << LinesUsed << " lines used for calibration ";
      for (i = 0; i < LinesUsed; i++) {
         cout << Energies[i] << " ";
      }
      cout << endl;
   }
   
   // Set up fit functions
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

   // And fit.
   Opts = FitOptions;
   if (LinesUsed > 1) {
      CalibPlot.Fit(CalibFitLin, Opts.c_str());
      HistoCal.LinGainFit[0] = CalibFitLin->GetParameter(0);
      HistoCal.dLinGainFit[0] = CalibFitLin->GetParError(0);
      HistoCal.LinGainFit[1] = CalibFitLin->GetParameter(1);
      HistoCal.dLinGainFit[1] = CalibFitLin->GetParError(1);
      HistoCal.LinGainFit[2] = CalibFitLin->GetChisquare() / CalibFitLin->GetNDF();
   }
   Opts += "+";
   if (LinesUsed > 2) {
      CalibPlot.Fit(CalibFitQuad, Opts.c_str());
      HistoCal.QuadGainFit[0] = CalibFitQuad->GetParameter(0);
      HistoCal.dQuadGainFit[0] = CalibFitQuad->GetParError(0);
      HistoCal.QuadGainFit[1] = CalibFitQuad->GetParameter(1);
      HistoCal.dQuadGainFit[1] = CalibFitQuad->GetParError(1);
      HistoCal.QuadGainFit[2] = CalibFitQuad->GetParameter(2);
      HistoCal.dQuadGainFit[2] = CalibFitQuad->GetParError(2);
      HistoCal.QuadGainFit[3] = CalibFitQuad->GetChisquare() / CalibFitQuad->GetNDF();
   }

   if (Config.PrintVerbose) {
      
      // Removed solution part from here because it is done at peak search stage.
      // Maybe pass the values in to here from fitting stage if this is actually needed.
      
      //cout << endl << "Two point solution:" << endl;
      //cout << "Offset = " << O << " +/- " << dO << "\tGain = " << G << " +/- " << dG << endl << endl;
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
   GainPlot->SetBinContent(ItemNum + 1, HistoCal.QuadGainFit[1]);  // ItemNum+1 to skip bin 0 which is underflow
   OffsetPlot->SetBinContent(ItemNum + 1, HistoCal.QuadGainFit[0]);
   QuadPlot->SetBinContent(ItemNum + 1, HistoCal.QuadGainFit[2]);

   GainHist->Fill(HistoCal.QuadGainFit[1]);
   OffsetHist->Fill(HistoCal.QuadGainFit[0]);
   QuadHist->Fill(HistoCal.QuadGainFit[2]);
   
   // Now print reports on results of fits and calibration.
   if (HistoCal.LinesUsed > 1) {
      if (HistoCal.LinesUsed < 3 || FORCE_LINEAR) {
         GainOut << Settings.OutputName << ":\t" << HistoCal.LinGainFit[0];
         GainOut << "\t" << HistoCal.LinGainFit[1] << endl;
      } else {
         GainOut << Settings.OutputName << ":\t" << HistoCal.QuadGainFit[0] << "\t";
         GainOut << HistoCal.QuadGainFit[1] << "\t" << HistoCal.QuadGainFit[2] << endl;
      }
   } else {
      //GainOut << HistName << " Fail!!!" << endl;
   }
   // Write full calibration report
   if (HistoCal.LinesUsed > 0) {
      CalibrationReport(&HistoCal, ReportOut, Settings.OutputName, Settings);
   } else {
      ReportOut << endl << "------------------------------------------" << endl << Settings.OutputName << endl <<
          "------------------------------------------" << endl << endl;
      ReportOut << "Fail Fail Fail! The calibration has failed!" << endl;
   }
   
   // Write .cal file for GRSISpoon
   if (Config.CalFile) {
      WriteCalFile(&HistoCal, CalFileOut, Settings.HistName, Settings);
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


int ConfigureEnergyFit(int Clover, int Crystal, int Seg,  int FileType, int FileNum, FitSettings *Settings) {
   
   char CharBuf[CHAR_BUFFER_SIZE];
   int x;
   
   // Set source and histogram name and output name based on seg number and file type
   // This should be factored out to a separate function
   if (FileType == 1) {
      switch (Seg) {
      case 0:
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da Chg", Clover, Num2Col(Crystal), Seg);
         Settings->HistName = CharBuf;
         Settings->OutputName = CharBuf;
         Settings->Source = Config.SourceNumCore.at(FileNum);
         break;
      case 9:
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db Chg", Clover, Num2Col(Crystal), 0);
         Settings->HistName = CharBuf;
         Settings->OutputName = CharBuf;
         Settings->Source = Config.SourceNumCore.at(FileNum);
         break;
      default:
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx Chg", Clover, Num2Col(Crystal), Seg);
         Settings->HistName = CharBuf;
         Settings->OutputName = CharBuf;
         if (Seg < 5) {
            Settings->Source = Config.SourceNumFront.at(FileNum);
         } else {
            Settings->Source = Config.SourceNumBack.at(FileNum);
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
         Settings->HistName = CharBuf;
         Settings->OutputName = CharBuf;
         Settings->OutputName += "\t";
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da Chg", Clover, Num2Col(Crystal), Seg);
         Settings->OutputName += CharBuf;
         Settings->Source = Config.SourceNumCore.at(FileNum);
         break;
      case 9:
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "Chrg0%03d", ((Clover - 1) * 60) + x + 9);
         Settings->HistName = CharBuf;
         Settings->OutputName = CharBuf;
         Settings->OutputName += "\t";
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db Chg", Clover, Num2Col(Crystal), 0);
         Settings->OutputName += CharBuf;
         Settings->Source = Config.SourceNumCore.at(FileNum);
         break;
      default:
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "Chrg0%03d", ((Clover - 1) * 60) + x + Seg);
         Settings->HistName = CharBuf;
         Settings->OutputName = CharBuf;
         Settings->OutputName += "\t";
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx Chg", Clover, Num2Col(Crystal), Seg);
         Settings->OutputName += CharBuf;
         if (Seg < 5) {
            Settings->Source = Config.SourceNumFront.at(FileNum);
         } else {
            Settings->Source = Config.SourceNumBack.at(FileNum);
         }
         break;
      }


   }
   
   // Check if plot should be active for this channel
   Settings->PlotOn = 0;
   if (Config.PlotFits) {
      if (Config.CalibPlots[Clover - 1][Crystal][Seg]) {
         Settings->PlotOn = 1;
      }
   }
   // Check if Manual peak select should be active
   if (Config.ManualPeakSelect[Clover - 1][Crystal][Seg]) {
      Settings->PeakSelect  = 1;
   }
   else {
      Settings->PeakSelect = 0;
   }
   
   if (FileType == 1) {
      Settings->Integration = INTEGRATION;
      Settings->Dispersion = float (CHARGE_BINS) / float (CHARGE_MAX);
   }
   if (FileType == 2) {  // histogram file from analyser
      Settings->Integration = 1;
      Settings->Dispersion = 1;
   }
   
   Settings->SearchSigma = EN_SEARCH_SIGMA;
   Settings->SearchThresh = EN_SEARCH_THRESH;
   Settings->SigmaEstZero = ENERGY_SIGMA_ZERO;
   Settings->SigmaEst1MeV = ENERGY_SIGMA_1MEV;
   Settings->FitZero = Config.FitZero;
   Settings->BackupPeakSelect = 0;
   
   return 0;
}

int ConfigureWaveEnFit(int Clover, int Crystal, int Seg,  int FileType, int FileNum, FitSettings *Settings) {
   
   char CharBuf[CHAR_BUFFER_SIZE];
   
   if (Seg == 0 || Seg == 9) {
      Settings->Source = Config.SourceNumCore.at(FileNum);
      if (Seg == 0) {
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da WaveChg", Clover, Num2Col(Crystal), Seg);
         Settings->HistName = CharBuf;
         Settings->OutputName = CharBuf;
      } else {
         snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db WaveChg", Clover, Num2Col(Crystal), 0);
         Settings->HistName = CharBuf;
         Settings->OutputName = CharBuf;
      }
   } else {
      snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx WaveChg", Clover, Num2Col(Crystal), Seg);
      Settings->HistName = CharBuf;
      Settings->OutputName = CharBuf;
      if (Seg < 5) {
         Settings->Source = Config.SourceNumFront.at(FileNum);
      } else {
         Settings->Source = Config.SourceNumBack.at(FileNum);
      }
   }
   
   // Check if plot should be active for this channel
   Settings->PlotOn = 0;
   if (Config.PlotFits) {
      if (Config.CalibPlots[Clover - 1][Crystal][Seg]) {
         Settings->PlotOn = 1;
      }
   }
   // Check if Manual peak select should be active
   if (Config.ManualPeakSelect[Clover - 1][Crystal][Seg]) {
      Settings->PeakSelect  = 1;
   }
   else {
      Settings->PeakSelect = 0;
   }
   
   Settings->Integration = 1;
   Settings->Dispersion = float (CHARGE_BINS) / float (WAVE_CHARGE_MAX);
   Settings->SearchSigma = WAVE_SEARCH_SIGMA;
   Settings->SearchThresh = WAVE_SEARCH_THRESH;
   Settings->SigmaEstZero = WAVE_SIGMA_ZERO;
   Settings->SigmaEst1MeV = WAVE_SIGMA_1MEV;
   Settings->FitZero = Config.FitZero;
   Settings->BackupPeakSelect = 0;
   
   return 0;
}

