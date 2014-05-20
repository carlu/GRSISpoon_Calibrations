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
#include <TSystem.h>
#include <TFolder.h>

// TriScope libraries
//#include "TTigFragment.h"

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "Calib.h"
#include "CalibTools.h"
#include "Utils.h"

extern TApplication *App;

// Canvas for plots
TCanvas *cCalib1, *cCalib1a, *cCalib2,  *ctemp;
// Output file stuff
TFile *outfile;
TDirectory *dCharge, *dWaveCharge, *dCalibration = { 0 };

// output streams 
ofstream GainOut;
ofstream ReportOut;
ofstream WaveOut;
ofstream WaveReportOut;
ofstream CalFileOut;

// Spectra for Calibration
TH1F *GainPlot, *OffsetPlot, *QuadPlot;
TH1F *GainHist, *OffsetHist, *QuadHist;

// Fitting stuff
std::string FitOptions = ("RQE");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace
std::string Opts;

// ------------------------------------------------
// Functions for managing calibration in SortHistos
// ------------------------------------------------

// Loop all files/sources
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
   
   // Check we have sources for all runs
   if(Config.files.size() != Config.SourceNumCore.size()) {
      cout << "For histogram calibration the number of files (" << Config.files.size() << ") must equal number of sources (" << Config.SourceNumCore.size()  << ")" << endl;
      return -1;
   }
      
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
      FitHistoFile(file, FileType, FileNum, &FitMap, &WaveFitMap);
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
               
               HistoFit Fit;
               HistoCal Cal;
               
               ConfigureEnergyFit(Clover,Crystal, Seg, FileType, FileNum, &Settings);
               
               if (Config.PrintVerbose) {
                  cout << endl << "------------------------------------------------------------------------" << endl;
                  cout << "Calibrating " << Settings.HistName << endl;
                  cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                  cout << "------------------------------------------------------------------------" << endl << endl;
               }
               
               CalibSuccess = CalibrateChannel(FitMap[ChanVector], Settings, &Fit, &Cal, GainOut, ReportOut);
            }
         }   
      }
   }
   
   if (Config.CalEnergy) {
      GainOut.close();
   }
   if (Config.CalReport) {
      ReportOut.close();
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
               
               HistoFit Fit;
               HistoCal Cal;
               
               ConfigureWaveEnFit(Clover,Crystal, Seg, FileType, FileNum, &Settings);
               
               if (Config.PrintVerbose) {
                  cout << endl << "------------------------------------------------------------------------" << endl;
                  cout << "Calibrating " << Settings.HistName << endl;
                  cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                  cout << "------------------------------------------------------------------------" << endl << endl;
               }
               
               CalibSuccess = CalibrateChannel(WaveFitMap[ChanVector], Settings, &Fit, &Cal, WaveOut, WaveReportOut);
            
            }
         }   
      }
   }
   
   if (Config.CalEnergy) {
      WaveOut.close();
   }
   if (Config.CalReport) {
      WaveReportOut.close();
   }
   
   if (Config.WriteFits) {
      outfile->Close();
   }
   
   return 0;
}

// Loop histograms in one file
int FitHistoFile(TFile *file, int FileType, int FileNum, MasterFitMap *FitMap, MasterFitMap *WaveFitMap) {

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
                     cout << endl << "------------------------------------------------------------------------" << endl;
                     cout << "Hist " << Settings.HistName << " loaded" << endl;
                     cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                     cout << "------------------------------------------------------------------------" << endl << endl;
                  }
                  
                  
                  
                  // Initialise FitSettings and HistoFit               
                  HistoFit Fit;
                  HistoCal Cal;
                  
                  // Perform fit
                  FitSuccess = FitGammaSpectrum(Histo, &Fit, &Cal, Settings);
                  
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
                  if(!FitMap->count(ChanVector)) {  // If entry for this channel does not exist...
                     FitMap->insert(MasterFitPair(ChanVector,ChanFits));  // create it....
                  }
                  else {                                                   // else....
                     
                     FitMap->at(ChanVector).insert(ChanFits.begin(),ChanFits.end());      // add to it.
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
               
               TH1F *Histo = (TH1F *) file->FindObjectAny(Settings.HistName.c_str());
               
               if(Histo) {
               
                  if (Config.PrintVerbose) {
                     cout << endl << "------------------------------------------------------------------------" << endl;
                     cout << "Hist " << Settings.HistName << " loaded" << endl;
                     cout << "Clover " << Clover << " Crystal " << Crystal << " Seg " << Seg << endl;
                     cout << "------------------------------------------------------------------------" << endl << endl;
                  }
                     
                  // Perform Fit                  
                  HistoFit WaveFit;
                  HistoCal WaveCal;
                  
                  FitSuccess = FitGammaSpectrum(Histo, &WaveFit, &WaveCal, Settings);
                  
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
                  if(!WaveFitMap->count(ChanVector)) {  // If entry for this channel does not exist...
                     WaveFitMap->insert(MasterFitPair(ChanVector,ChanFits));  // create it....
                  }
                  else {                                                   // else....
                     
                     WaveFitMap->at(ChanVector).insert(ChanFits.begin(),ChanFits.end());      // add to it.
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

   return 0;
}

// --------------------------------
// Functions for fitting a spectrum 
// --------------------------------

// Find/id first two peaks, rough calibration, loop all peaks
int FitGammaSpectrum(TH1F * Histo, HistoFit * Fit, HistoCal * Cal, FitSettings Settings)
{

   // Histo:         pointer to histo to be fitted
   // Fit:           HistoFit structure into which the fit results will be written
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
   int Line; //, LinesUsed;
   float G, O, dG, dO, FitCentre;       // values for fast 2 point calibration.
   float En1, En2;              // line energies for two point calib
   int CustomPeak1;
   int CustomPeak2;
   float CustomCentroid1;
   float CustomCentroid2;
   En1 = Config.Sources[Settings.Source][0];
   En2 = Config.Sources[Settings.Source][1];
   float IdealRatio = En1 / En2;
   float Chg1, Chg2, dChg1, dChg2;      // charge and erros for 2 point
   FitResult FitRes[MAX_LINES+1]; // store full fit results
   int Integral;                // counts in spectrum
   TSpectrum *Spec = new TSpectrum();
   TF1 *FitRange[MAX_LINES];    // pointers to functions for fitting each line
   float InitialGain = 0.16;
   // Fitting stuff
   std::string FitOptions = ("RQEM");
   // R=restrict to function range, Q=quiet, L=log likelihood method, E=improved err estimation, + add fit instead of replace

   if (Histo) {
      Integral = Histo->Integral();
   } else {
      cout << "Bad histo! " << endl;
   }

   if (Config.PrintVerbose) {
      cout << "\tIntegral: " << Integral;
   }
   if (Integral > MIN_FIT_COUNTS) {

      if (Settings.PlotOn || Settings.PeakSelect || Settings.BackupPeakSelect) {
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
                  Diff = fabs(Ratio - IdealRatio) / IdealRatio;
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
            cout << "BestPeak1: " << BestPeak1 << " BestPeak2: " << BestPeak2 << " Ratio: " << PeakPositions[BestPeak1]
                / PeakPositions[BestPeak2] << endl;
         }
      } else {
         return -2;
      }
       
      // Plot spectrum 
      if(Settings.PlotOn || Settings.PeakSelect || (PeakFound==0 && Settings.BackupPeakSelect)) {
         cCalib1->cd(1);
         Histo->Draw();
         //cCalib1->Modified();
         cCalib1->Update();
         App->Run(1);
      }
      
      // Manual Peak Selection
      if (Settings.PeakSelect || (PeakFound==0 && Settings.BackupPeakSelect)) {
         cout << "Peaks: " << endl;
         for(i=0; i<NumPeaks; i++) {
            cout << i << ":\t" << PeakPositions[i];
            if(i==BestPeak1) {
               cout << "  **Peak1 (" << En1 << " keV)";
            }
            if(i==BestPeak2) {
               cout << "  **Peak2 (" << En2 << " keV)";
            }
            cout << endl;
         }
         
         //gSystem->ProcessEvents();
         Histo->Draw();
         cCalib1->Update();
         App->Run(1);

         cout << "Please enter peak number for " << En1 << " keV (-1 to specify custom centroid)" << endl;
         cin >> CustomPeak1;
         if(CustomPeak1>NumPeaks) {
            cout << "What?!" << endl;
            return -3;
         }
         if(CustomPeak1<0){
            cout << "Enter custom centroid:" << endl;
            cin >> CustomCentroid1;
         }
         else {
            CustomCentroid1 = PeakPositions[CustomPeak1];
         }
         
         cout << "Please enter peak number for " << En2 << " keV (-1 to specify custom centroid)" << endl;
         cin >> CustomPeak2;
         if(CustomPeak2>NumPeaks) {
            cout << "What?!" << endl;
            return -3;
         }
         if(CustomPeak2<0){
            cout << "Enter custom centroid:" << endl;
            cin >> CustomCentroid2;
         }
         else {
            CustomCentroid2 = PeakPositions[CustomPeak2];
         }
         
         PeakPositions[0] = CustomCentroid1;
         BestPeak1 = 0;
         PeakPositions[1] = CustomCentroid2;
         BestPeak2 = 1;
         PeakFound = 1;
         cout << "Proceeding with calibration using the following peaks: " << endl;
         cout << "\t" << PeakPositions[0] << " ch = " << En1 << " keV" << endl;
         cout << "\t" << PeakPositions[1] << " ch = " << En2 << " keV" << endl;
      
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
            
            // Store fit result for use in calibration                      
            //memcpy(&Fit->PeakFits[Line], &FitRes[Line], sizeof(FitResult));
            Fit->PeakFits.push_back(FitRes[Line]);
            Fit->FitSuccess.push_back(1);
         }

         //-------------------------------------------------------------//
         // Calculate a rough calibration                               //
         //-------------------------------------------------------------//

         Chg1 = FitRes[0].Mean / Settings.Integration;
         Chg2 = FitRes[1].Mean / Settings.Integration;
         dChg1 = FitRes[0].dMean / Settings.Integration;
         dChg2 = FitRes[1].dMean / Settings.Integration;

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

         Cal->LinGain[0] = O;
         Cal->dLinGain[0] = dO;
         Cal->LinGain[1] = G;
         Cal->dLinGain[1] = dG;

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
               
               // Store fit result for use in calibration                      
               //memcpy(&Fit->PeakFits[Line], &FitRes[Line], sizeof(FitResult));
               //Fit->FitSuccess[Line] = 1;
               Fit->PeakFits.push_back(FitRes[Line]);
               Fit->FitSuccess.push_back(1);
            }
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

   return 0;

}

// Fit a peak
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


   if (Config.PrintVerbose) {
      cout << "Fitting line " << Line << " Min: " << Min << " Max: " << Max << endl;
   }
   if (Settings.PlotOn) {
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

// --------------------------
// Functions for calibration
// --------------------------

// Calibrate a channel
int CalibrateChannel(ChannelFitMap Fits, FitSettings Settings, HistoFit *Fit, HistoCal * Cal, ofstream &Out, ofstream &Report ) {
   
   // This function should copy the operation of the old CalibrateGammaSpectrum() in CalibTools.C
   // However it should be altered to use maps and vectors rather than arrays.
   int Clover, Crystal, Seg, i;
   int LinesUsed;
   int ItemNum;
   bool TestsPassed;
   std::string Name;
   float Energy = 0.0;
   // Would like the following to be vectors but can't get the TGraphErrors to work 
   float Energies[MAX_TOTAL_LINES], dEnergies[MAX_TOTAL_LINES];
   float Centroids[MAX_TOTAL_LINES], dCentroids[MAX_TOTAL_LINES];
   
   FitResult LineFit;
   
   LinesUsed = 0;
   
   for (ChannelFitMapIt Line = Fits.begin(); Line != Fits.end(); Line++) {
      TestsPassed = 0;
      Energy = Line->first;
      LineFit = Line->second;
      
      Fit->PeakFits.push_back(LineFit); // Store fit for use in calibration and report 
      
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
            Fit->FitSuccess.push_back(1);
            if(LinesUsed > (MAX_TOTAL_LINES -2)) {  // Stop adding lines if we are going to overflow arrays
               break;                                // -1 here, -2 because of poss of 0ch=0keV line
            }
            
         }
      }
      if(TestsPassed==0) {
         Fit->FitSuccess.push_back(0);
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
      Fit->PeakFits.push_back(ZeroFit);
      Fit->FitSuccess.push_back(1);
   }   
   
   Cal->LinesUsed = LinesUsed;
   
   // Now add points to a TGraph 
   TGraphErrors CalibPlot(LinesUsed, Centroids, Energies, dCentroids, dEnergies);
   Name = Settings.HistName + " Calibration";
   CalibPlot.SetTitle(Name.c_str());
   Name = Settings.HistName + " Cal";
   CalibPlot.SetName(Name.c_str());
   CalibPlot.GetXaxis()->SetTitle("Centroid (ch)");
   CalibPlot.GetYaxis()->SetTitle("Energy (keV)");
   CalibPlot.GetYaxis()->SetLimits(0.0,2000.0);
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
      Cal->LinGainFit[0] = CalibFitLin->GetParameter(0);
      Cal->dLinGainFit[0] = CalibFitLin->GetParError(0);
      Cal->LinGainFit[1] = CalibFitLin->GetParameter(1);
      Cal->dLinGainFit[1] = CalibFitLin->GetParError(1);
      Cal->LinGainFit[2] = CalibFitLin->GetChisquare() / CalibFitLin->GetNDF();
   }
   Opts += "+";
   if (LinesUsed > 2) {
      CalibPlot.Fit(CalibFitQuad, Opts.c_str());
      Cal->QuadGainFit[0] = CalibFitQuad->GetParameter(0);
      Cal->dQuadGainFit[0] = CalibFitQuad->GetParError(0);
      Cal->QuadGainFit[1] = CalibFitQuad->GetParameter(1);
      Cal->dQuadGainFit[1] = CalibFitQuad->GetParError(1);
      Cal->QuadGainFit[2] = CalibFitQuad->GetParameter(2);
      Cal->dQuadGainFit[2] = CalibFitQuad->GetParError(2);
      Cal->QuadGainFit[3] = CalibFitQuad->GetChisquare() / CalibFitQuad->GetNDF();
   }
   
   if(Config.WriteFits && !Settings.TempFit) {
      dCalibration->cd();
      CalibPlot.Write();
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


   if (Config.PlotCalib && Settings.PlotOn && !Settings.TempFit) {
      cCalib1a->cd(1);
      CalibPlot.SetMarkerColor(2);
      CalibPlot.SetMarkerStyle(20);
      CalibPlot.SetMarkerSize(1.0);
      CalibPlot.SetTitle("Calibration");
      CalibPlot.Draw("AP");
      App->Run(1);
      //App->Run();
   }
   
   // generate output....
   
   if(!Settings.TempFit) {
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
      
      
      GainPlot->SetBinContent(ItemNum + 1, Cal->QuadGainFit[1]);  // ItemNum+1 to skip bin 0 which is underflow
      OffsetPlot->SetBinContent(ItemNum + 1, Cal->QuadGainFit[0]);
      QuadPlot->SetBinContent(ItemNum + 1, Cal->QuadGainFit[2]);

      GainHist->Fill(Cal->QuadGainFit[1]);
      OffsetHist->Fill(Cal->QuadGainFit[0]);
      QuadHist->Fill(Cal->QuadGainFit[2]);
   }
   
   // Now print reports on results of fits and calibration.
   if (Cal->LinesUsed > 1) {
      if (Cal->LinesUsed < 3 || FORCE_LINEAR) {
         Out << Settings.OutputName << ":\t" << Cal->LinGainFit[0];
         Out << "\t" << Cal->LinGainFit[1] << endl;
      } else {
         Out << Settings.OutputName << ":\t" << Cal->QuadGainFit[0] << "\t";
         Out << Cal->QuadGainFit[1] << "\t" << Cal->QuadGainFit[2] << endl;
      }
   } else {
      //GainOut << HistName << " Fail!!!" << endl;
   }
   // Write full calibration report
   if (Cal->LinesUsed > 0) {
      CalibrationReport(Fit, Cal, Report, Settings.OutputName, Settings);
   } else {
      Report << endl << "------------------------------------------" << endl << Settings.OutputName << endl <<
          "------------------------------------------" << endl << endl;
      Report << "Fail Fail Fail! The calibration has failed!" << endl;
   }
   
   // Write .cal file for GRSISpoon
   if (Config.CalFile) {
      WriteCalFile(Fit, CalFileOut, Settings.HistName, Settings);
   }

   


   if (Config.PlotCalibSummary) {
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
   
   
   return 0;

}

// Other general helper functions
// ------------------------------
// Configure fit settings for Energy/WaveEn spectrum
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
   Settings->TempFit = 0;
   
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
   Settings->TempFit = 0;
   
   return 0;
}

// Generate text report of fits and calibrations
int CalibrationReport(HistoFit * Fit, HistoCal * Cal, ofstream & ReportOut, std::string HistName, FitSettings Settings)
{

   int i, NumFits;
   float Energies[MAX_LINES], Residuals[MAX_LINES];
   float CalibEn;

   if(Config.CalReport) {
      // Heading for this channel
      ReportOut << endl << "------------------------------------------" << endl << HistName << endl <<
          "------------------------------------------" << endl << endl;
      // Peak fit info
      ReportOut << Cal->LinesUsed << " lines used in calibration (*)" << endl;
      ReportOut << "Individual fit results:" << endl;
      ReportOut << "(Const +/- err\tMean +/- err\tSigma +/- err\nChiSq\tNDF\tCSPD)" << endl << endl;
      //for (i = 0; i <= Config.Sources[Settings.Source].size(); i++) {
      for (i = 0; i < Fit->PeakFits.size(); i++) {  
         ReportOut << Fit->PeakFits.at(i).Energy << " keV ";
         if (Fit->FitSuccess.at(i) == 1) {
            ReportOut << "(*)";
         }
         ReportOut << endl;
         ReportOut << Fit->PeakFits.at(i).Const << " +/- " << Fit->PeakFits.at(i).dConst << "\t";
         ReportOut << Fit->PeakFits.at(i).Mean << " +/- " << Fit->PeakFits.at(i).dMean << "\t";
         ReportOut << Fit->PeakFits.at(i).Const << " +/- " << Fit->PeakFits.at(i).dConst << endl;
         ReportOut << Fit->PeakFits.at(i).ChiSq << "\t" << Fit->PeakFits.at(i).NDF << "\t" << Fit->PeakFits.at(i).ChiSq /
             Fit->PeakFits.at(i).NDF;
         ReportOut << endl << endl;
      }
      // Calibration...
      ReportOut << "Linear Solution: Offset = " << Cal->LinGain[0] << " +/- " << Cal->dLinGain[0] << "\t";
      ReportOut << "Gain = " << Cal->LinGain[1] << " +/- " << Cal->dLinGain[1] << endl;
      ReportOut << "Linear Fit: Offset = " << Cal->LinGainFit[0] << " +/- " << Cal->dLinGainFit[0] << "\t";
      ReportOut << "Gain = " << Cal->LinGainFit[1] << " +/- " << Cal->dLinGainFit[1] << "\t";
      ReportOut << "CSPD = " << Cal->LinGainFit[2] << endl;
      ReportOut << "Quadratic Fit: Offset = " << Cal->QuadGainFit[0] << " +/- " << Cal->dQuadGainFit[0] << "\t";
      ReportOut << "Gain = " << Cal->QuadGainFit[1] << " +/- " << Cal->dQuadGainFit[1] << "\t";
      ReportOut << "Quad = " << Cal->QuadGainFit[2] << " +/- " << Cal->dQuadGainFit[2] << "\t";
      ReportOut << "CSPD = " << Cal->LinGainFit[3] << endl;
      // Residual from quadratic fit
      ReportOut << endl << "Quadratic calibration residuals...." << endl;
      ReportOut << "Centroid (ch)\t\tList Energy (keV)\t\tCalibration Energy (keV)\t\tResidual (keV)" << endl;
      NumFits = 0;
      //for (i = 0; i <= Config.Sources[Settings.Source].size(); i++) {
      for (i = 0; i < Fit->PeakFits.size(); i++) {
         if (Fit->FitSuccess.at(i) == 1) {
            CalibEn = Cal->QuadGainFit[0] + (Cal->QuadGainFit[1] * (Fit->PeakFits.at(i).Mean / Settings.Integration)) +
                (pow((Fit->PeakFits.at(i).Mean / Settings.Integration), 2) * Cal->QuadGainFit[2]);
            ReportOut << Fit->PeakFits.at(i).Mean << "\t\t\t" << Fit->PeakFits.at(i).Energy << "\t\t\t";
            ReportOut << CalibEn << "\t\t\t" << CalibEn - Fit->PeakFits.at(i).Energy << endl;
            Residuals[NumFits] = CalibEn - Fit->PeakFits.at(i).Energy;
            NumFits++;
         }
      }
      // Residual from linear fit
      NumFits = 0;
      ReportOut << endl << "Linear calibration residuals...." << endl;
      ReportOut << "Centroid (ch)\t\tList Energy (keV)\t\tCalibration Energy (keV)\t\tResidual (keV)" << endl;
      for (i = 0; i < Fit->PeakFits.size(); i++) {
      //for (i = 0; i <= Config.Sources[Settings.Source].size(); i++) {
         if (Fit->FitSuccess.at(i) == 1) {
            CalibEn = Cal->LinGainFit[0] + (Cal->LinGainFit[1] * (Fit->PeakFits.at(i).Mean / Settings.Integration));
            ReportOut << Fit->PeakFits.at(i).Mean << "\t\t\t" << Fit->PeakFits.at(i).Energy << "\t\t\t";
            ReportOut << CalibEn << "\t\t\t" << CalibEn - Fit->PeakFits.at(i).Energy << endl;
            Energies[NumFits] = Fit->PeakFits.at(i).Energy;
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

// Write out .cal file (not implemented yet)
int WriteCalFile(HistoFit *Fit, ofstream &CalFileOut, std::string HistName, FitSettings Settings) {



   return 0;
}
