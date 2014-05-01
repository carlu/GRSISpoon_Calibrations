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

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "Calib.h"
#include "CalibTools.h"

// File pointers:
static TFile *outfile = 0;
static TDirectory *dCharge, *dWaveCharge, *dTemp, *dOther = { 0 };

extern TApplication *App;

// Spectra Pointers:
static TH1F *hCharge[CLOVERS][CRYSTALS][SEGS + 2] = { };        // charge from FPGA
static TH1F *hWaveCharge[CLOVERS][CRYSTALS][SEGS + 2] = { };    // charge from waveform

static TH1F *hMidasTime = 0;
static TH1F *hCrystalChargeTemp[CLOVERS][CRYSTALS] = { };       // Only doing these guys for the cores
static TH1F *hCrystalGain[CLOVERS][CRYSTALS] = { };
static TH1F *hCrystalOffset[CLOVERS][CRYSTALS] = { };

static TH1F *WaveHist = 0;

// Functions called from main:   
void InitCalib();
void Calib(std::vector < TTigFragment > &ev);
void FinalCalib();
// Functions called from here:
void ResetTempSpectra();

extern TCanvas *cCalib1, *cCalib1a, *cCalib2, *cWave1, *ctemp;
// Storing run settings
extern RunConfig Config;

void InitCalib()
{

   char Colours[] = "BGRW";

   // Initialise output file                
   std::string tempstring = Config.OutPath + Config.CalOut;
   outfile = new TFile(tempstring.c_str(), "RECREATE");

   dCharge = outfile->mkdir("Charge");
   dWaveCharge = outfile->mkdir("WaveCharge");
   dTemp = outfile->mkdir("Temp");
   dOther = outfile->mkdir("Other");

   char name[512], title[512];
   int Clover, Crystal, Seg, NumBins;


   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         // temp charge histo
         dTemp->cd();
         sprintf(name, "TIG%02d%c Tmp Chg", Clover, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Temp Core Charge (arb)", Clover, Colours[Crystal]);
         hCrystalChargeTemp[Clover - 1][Crystal] = new TH1F(name, title, CHARGE_BINS, 0, CHARGE_MAX);
         // histo for record of calibration from temp spectra
         dOther->cd();
         sprintf(name, "TIG%02d%c Gain", Clover, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Gain vs Time", Clover, Colours[Crystal]);
         hCrystalGain[Clover - 1][Crystal] = new TH1F(name, title, TIME_BINS, 0, MAX_TIME);
         sprintf(name, "TIG%02d%c Offset", Clover, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Offset vs Time", Clover + 1, Colours[Crystal]);
         hCrystalOffset[Clover - 1][Crystal] = new TH1F(name, title, TIME_BINS, 0, MAX_TIME);
         // Spectra for charge from FPGA
         dCharge->cd();
         Seg = 0;
         sprintf(name, "TIG%02d%cN%02da Chg", Clover, Colours[Crystal], Seg);
         sprintf(title, "TIG%02d%cN%02da Charge (arb)", Clover, Colours[Crystal], Seg);
         hCharge[Clover - 1][Crystal][Seg] = new TH1F(name, title, CHARGE_BINS, 0, CHARGE_MAX);
         sprintf(name, "TIG%02d%cN%02db Chg", Clover, Colours[Crystal], Seg);
         sprintf(title, "TIG%02d%cN%02db Charge (arb)", Clover, Colours[Crystal], Seg);
         hCharge[Clover - 1][Crystal][9] = new TH1F(name, title, CHARGE_BINS, 0, CHARGE_MAX);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%cP%02dx Chg", Clover, Colours[Crystal], Seg);
            sprintf(title, "TIG%02d%cP%02dx Charge (arb)", Clover, Colours[Crystal], Seg);
            hCharge[Clover - 1][Crystal][Seg] = new TH1F(name, title, CHARGE_BINS, 0, CHARGE_MAX);
         }
         // and charge derived from waveform
         dWaveCharge->cd();
         Seg = 0;
         sprintf(name, "TIG%02d%cN%02da WaveChg", Clover, Colours[Crystal], Seg);
         sprintf(title, "TIG%02d%cN%02da Waveform Charge (arb)", Clover, Colours[Crystal], Seg);
         hWaveCharge[Clover - 1][Crystal][Seg] = new TH1F(name, title, CHARGE_BINS, 0, WAVE_CHARGE_MAX);
         sprintf(name, "TIG%02d%cN%02db WaveChg", Clover, Colours[Crystal], Seg);
         sprintf(title, "TIG%02d%cN%02db Waveform Charge (arb)", Clover, Colours[Crystal], Seg);
         hWaveCharge[Clover - 1][Crystal][9] = new TH1F(name, title, CHARGE_BINS, 0, WAVE_CHARGE_MAX);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%cP%02dx WaveChg", Clover, Colours[Crystal], Seg);
            sprintf(title, "TIG%02d%cP%02dx Waveform Charge (arb)", Clover, Colours[Crystal], Seg);
            hWaveCharge[Clover - 1][Crystal][Seg] = new TH1F(name, title, CHARGE_BINS, 0, WAVE_CHARGE_MAX);
         }
      }
   }
   // Time stamp histo
   dOther->cd();
   sprintf(name, "Midas Time");
   sprintf(title, "Midas Timestamps (s)");
   hMidasTime = new TH1F(name, title, TIME_BINS, 0, MAX_TIME);

   if (PLOT_WAVE) {
      sprintf(name, "Wavetemp");
      sprintf(title, "Temporary wave histogram");
      WaveHist = new TH1F(name, title, Config.WaveformSamples, 0, Config.WaveformSamples);
   }

   if (Config.PrintVerbose) {
      cout << "Searching for core Peaks: " << Config.Sources[Config.SourceNumCore][0] << "kev and " << Config.
          Sources[Config.SourceNumCore]
          [1] << "keV (Ratio " << Config.Sources[Config.SourceNumCore][0] /
          Config.Sources[Config.SourceNumCore][1] << ")" << endl;
      cout << "Searching for front seg peaks: " << Config.Sources[Config.SourceNumFront][0] << "kev and " << Config.
          Sources[Config.SourceNumFront][1] << "keV (Ratio " << Config.Sources[Config.SourceNumFront][0] /
          Config.Sources[Config.SourceNumFront][1] << ")" << endl;
      cout << "Searching for back seg peaks: " << Config.Sources[Config.SourceNumBack][0] << "kev and " << Config.
          Sources[Config.SourceNumBack][1] << "keV (Ratio " << Config.Sources[Config.SourceNumBack][0] /
          Config.Sources[Config.SourceNumBack][1] << ")" << endl;
   }

   if (PLOT_WAVE) {
      cWave1 = new TCanvas();
   }

}


void Calib(std::vector < TTigFragment > &ev)
{

   //Variables
   int i, j, k;
   int Samp, Length;
   int Slave, Port, Chan;
   int Crystal, Clover, SpecNum;
   int FitSuccess = 0;
   std::string name;

   time_t MidasTime;
   static time_t StartTime;
   static int FirstEvent = 1;
   double RunTimeElapsed;
   static double FitTimeElapsed = 0.0;

   int TimeBin = 0;
   float TB = 0.0;

   int PlotOn = 0;
   int Source = 0;

   float WaveCharge = 0.0;
   float WaveEnergy = 0.0;


   if (DEBUG) {
      cout << "--------- New Event ---------" << endl;
   }

   for (i = 0; i < ev.size(); i++) {

      //Get time of first fragment
      if (FirstEvent == 1) {
         StartTime = ev[i].MidasTimeStamp;      //ev[i].MidasTimeStamp;
         FirstEvent = 0;
         if (Config.PrintBasic) {
            cout << "MIDAS time of first fragment: " << ctime(&StartTime) << endl;
         }
      }
      // Insert test for time earlier than StartTime.....
      if (ev[i].MidasTimeStamp < StartTime) {
         StartTime = ev[i].MidasTimeStamp;      //ev[i].MidasTimeStamp;
         if (Config.PrintBasic) {
            cout << "Earlier event found! Updating time of first fragment to: " << ctime(&StartTime) << endl;
         }
      }

      Slave = ((ev[i].ChannelAddress & 0x00F00000) >> 20);
      Port = ((ev[i].ChannelAddress & 0x00000F00) >> 8);
      Chan = (ev[i].ChannelAddress & 0x000000FF);
      name = ev[i].ChannelName;
      //cout << "Slave, Port, Chan = " << Slave << ", " << Port << ", " << Chan << "\t" << name << endl;

      Mnemonic mnemonic;
      if (name.size() >= 10) {
         ParseMnemonic(&name, &mnemonic);
      } else {
         if (Config.PrintVerbose) {
            cout << "Bad mnemonic size: This shouldn't happen if the odb is correctly configured!" << endl;
         }
         continue;
      }

      // If TIGRESS HPGe then fill energy spectra
      if (mnemonic.system == "TI" && mnemonic.subsystem == "G") {
         // Determine Crystal
         char Colour = mnemonic.arraysubposition.c_str()[0];
         Crystal = Col2Num(Colour);
         if (Crystal == -1) {
            if (Config.PrintBasic) {
               cout << "Bad Colour: " << Colour << endl;
            }
            continue;
         }
         // Determine Clover position
         Clover = mnemonic.arrayposition;

         // Calcualte wave energy
         Length = ev[i].wavebuffer.size();
         if (Length > Config.WaveInitialSamples + Config.WaveFinalSamples) {
            //cout << name;// << endl;
            //cout << " samples: " << Length << endl;  
            WaveCharge = CalcWaveCharge(ev[i].wavebuffer);
            //cout << "Charge: " << WaveCharge << endl;
            if (PLOT_WAVE) {
               cWave1->cd();
               for (Samp = 0; Samp < Length; Samp++) {
                  WaveHist->SetBinContent(Samp, ev[i].wavebuffer.at(Samp));
               }

               WaveHist->Draw();
               cWave1->Modified();
               cWave1->Update();
               App->Run(1);
               //delete Wavetemp;
            }
         }
         // If Core
         if (mnemonic.segment == 0) {
            //cout << mnemonic.outputsensor << endl;
            if (Chan == 0) {
               if (!(mnemonic.outputsensor.compare("a") == 0 || mnemonic.outputsensor.compare("A") == 0)) {     // if this is the primary core output
                  if (Config.PrintBasic) {
                     cout << "Core Channel/Label mismatch" << endl;
                  }
               }
               if (ev[i].Charge > 0) {
                  // Increment raw charge spectra
                  if (DEBUG) {
                     cout << "A: Filling " << Clover
                         << ", " << Crystal << ", 0, " << mnemonic.outputsensor << " with charge = " << ev[i].
                         Charge << endl;
                  }
                  hCharge[Clover - 1][Crystal][0]->Fill(ev[i].Charge);
                  hWaveCharge[Clover - 1][Crystal][0]->Fill(WaveCharge);
               }
            } else {
               if (Chan == 9) {
                  if (!(mnemonic.outputsensor.compare("b") == 0 || mnemonic.outputsensor.compare("B") == 0)) {  // if this is the primary core output
                     if (Config.PrintBasic) {
                        cout << "Core Channel/Label mismatch" << endl;
                     }
                  }
                  if (ev[i].Charge > 0) {
                     // Increment raw charge spectra
                     if (DEBUG) {
                        cout << "B: Filling " << Clover << ", " << Crystal << ", 0, " << mnemonic.outputsensor <<
                            " with charge = " << ev[i].Charge << endl;
                     }
                     hCharge[Clover - 1][Crystal][9]->Fill(ev[i].Charge);
                     hWaveCharge[Clover - 1][Crystal][9]->Fill(WaveCharge);
                  }
               }
            }
         } else {
            if (mnemonic.segment < 9) {
               if (ev[i].Charge > 0) {
                  hCharge[Clover - 1][Crystal][mnemonic.segment]->Fill(ev[i].Charge);   // Fill segment spectra
                  hWaveCharge[Clover - 1][Crystal][mnemonic.segment]->Fill(WaveCharge);
               }
            }
         }
         // Now Get time elapsed in run
         MidasTime = ev[i].MidasTimeStamp;
         //cout << "Time: " << ctime(&MidasTime) << endl;
         RunTimeElapsed = difftime(MidasTime, StartTime);
         hMidasTime->Fill(RunTimeElapsed);
         //cout << "Time: " << ctime(&MidasTime) << endl;
         // Have we moved on to a new time period?
         if ((RunTimeElapsed - FitTimeElapsed >= TIME_BIN_SIZE) && FIT_TEMP_SPECTRA) {

            TB = ((MidasTime - StartTime) / MAX_TIME);
            TimeBin = TB * TIME_BINS;

            if (DEBUG) {
               cout << "MT: " << MidasTime << " ST: " << StartTime << " MT: " << MAX_TIME << " TIME_BINS: " << TIME_BINS
                   << " TB: " << TB << " TimeBin: " << TimeBin << endl;
            }

            if (Config.PrintVerbose) {
               cout << "Start Of Run:      " << ctime(&StartTime);
               cout << "Start of this fit: " << FitTimeElapsed << endl;
               cout << "Time of this frag: " << ctime(&MidasTime);
               cout << "RunTimeElapsed:" << RunTimeElapsed << endl;
               cout << "TimeBin: " << TimeBin << endl;
            }
            FitTimeElapsed = RunTimeElapsed;


            // Loop crystals, fit spectra core , write gains and offsets to spectrum
            for (j = 1; j <= CLOVERS; j++) {
               for (k = 0; k < CRYSTALS; k++) {
                  if (Config.PrintVerbose) {
                     cout << "Clov: " << j << " Crys: " << k;
                  }
                  // Perform Fit                  
                  SpectrumFit Fit = { 0 };
                  FitSettings Settings = { 0 };

                  Settings.Source = Config.SourceNumCore;
                  Settings.Integration = INTEGRATION;
                  Settings.Dispersion = float (CHARGE_BINS) / float (CHARGE_MAX);
                  Settings.SearchSigma = EN_SEARCH_SIGMA;
                  Settings.SearchThresh = EN_SEARCH_THRESH;
                  Settings.SigmaEstZero = ENERGY_SIGMA_ZERO;
                  Settings.SigmaEst1MeV = ENERGY_SIGMA_1MEV;
                  Settings.FitZero = INCLUDE_ZERO;
                  Settings.PlotOn = 0;  // don't plot temp spectra fits
                  Settings.PeakSelect = 0;
                  Settings.BackupPeakSelect = 0;
                  
                  FitSuccess = FitGammaSpectrum(hCharge[Clover - 1][Crystal][0], &Fit, Settings);

                  if (FitSuccess > -1) {
                     hCrystalGain[j - 1][k]->SetBinContent(TimeBin, Fit.LinGainFit[1]);
                     hCrystalOffset[j - 1][k]->SetBinContent(TimeBin, Fit.LinGainFit[0]);
                  } else {
                     if (VERBOSE) {
                        switch (FitSuccess) {
                        case -1:
                           continue;
                        case -2:
                           continue;
                        case -3:
                           continue;
                        case -4:
                           continue;
                        default:
                           continue;
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
      // If TIGRESS suppressor
      if (mnemonic.system == "TI" && mnemonic.subsystem == "G") {

      }
   }
}



void FinalCalib()
{

   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   int FitSuccess = 0;
   bool PlotOn = 0;
   bool PeakSelect=0;
   int Source = 0;
   int i;
   ofstream GainOut;
   ofstream ReportOut;
   ofstream CalFileOut;
   ofstream WaveOut;
   ofstream WaveReportOut;
   string HistName;
   char CharBuf[CHAR_BUFFER_SIZE];

   TH1F *GainPlot, *OffsetPlot, *QuadPlot;
   TH1F *GainHist, *OffsetHist, *QuadHist;
   int ItemNum;


   float CalibEn = 0.0;

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

   std::string tempstring;
   // Prepare energy calibration output files
   if (Config.CalEnergy) {
      tempstring = Config.OutPath + "GainsOut.txt";
      GainOut.open(tempstring.c_str());
      if (Config.CalReport) {
         tempstring = Config.OutPath + "CalibrationReport.txt";
         ReportOut.open(tempstring.c_str());
      }
   }
   // Prepare wave calibration output files
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
   // Write spectra to file
   outfile->cd();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg <= SEGS + 1; Seg++) {
            dCharge->cd();
            hCharge[Clover - 1][Crystal][Seg]->Write();
            dWaveCharge->cd();
            hWaveCharge[Clover - 1][Crystal][Seg]->Write();
         }
         dTemp->cd();
         hCrystalChargeTemp[Clover - 1][Crystal]->Write();
         dOther->cd();
         hCrystalGain[Clover - 1][Crystal]->Write();
         hCrystalOffset[Clover - 1][Crystal]->Write();
      }
   }
   hMidasTime->Write();

   // Fit full-run energy spectra
   // -----------------------------
   if (Config.CalEnergy) {
      if (Config.PrintVerbose) {
         cout << endl << "Now fitting energy spectra from whole run..." << endl;
      }

      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            for (Seg = 0; Seg <= SEGS + 1; Seg++) {

               // Set source and histogram name based on seg number
               switch (Seg) {
               case 0:
                  snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da Chg", Clover, Num2Col(Crystal), Seg);
                  HistName = CharBuf;
                  Source = Config.SourceNumCore;
                  break;
               case 9:
                  snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db Chg", Clover, Num2Col(Crystal), 0);
                  HistName = CharBuf;
                  Source = Config.SourceNumCore;
                  break;
               default:
                  snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx Chg", Clover, Num2Col(Crystal), Seg);
                  HistName = CharBuf;
                  if (Seg < 5) {
                     Source = Config.SourceNumFront;
                  } else {
                     Source = Config.SourceNumBack;
                  }
                  break;
               }

               // Load histogram
               TH1F *Histo = hCharge[Clover - 1][Crystal][Seg];
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
                  // Perform Fit                  
                  SpectrumFit Fit = { 0 };
                  FitSettings Settings = { 0 };

                  Settings.Source = Source;
                  Settings.Integration = INTEGRATION;
                  Settings.Dispersion = float (CHARGE_BINS) / float (CHARGE_MAX);
                  Settings.SearchSigma = EN_SEARCH_SIGMA;
                  Settings.SearchThresh = EN_SEARCH_THRESH;
                  Settings.SigmaEstZero = ENERGY_SIGMA_ZERO;
                  Settings.SigmaEst1MeV = ENERGY_SIGMA_1MEV;
                  Settings.FitZero = Config.FitZero;
                  Settings.PlotOn = PlotOn;
                  Settings.PeakSelect = PeakSelect;
                  Settings.BackupPeakSelect = 0;
                  
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
                  if (OUTPUT_GAIN) {
                     if (FitSuccess > 0) {
                        if (FitSuccess < 3 || FORCE_LINEAR) {
                           GainOut << HistName << ":\t" << Fit.LinGainFit[0];
                           GainOut << "\t" << Fit.LinGainFit[1] << endl;
                        } else {
                           GainOut << HistName << ":\t" << Fit.QuadGainFit[0] << "\t";
                           GainOut << Fit.QuadGainFit[1] << "\t" << Fit.QuadGainFit[2] << endl;
                        }
                     } else {
                        //GainOut << HistName << " Fail!!!" << endl;
                     }
                  }
                  // Write full calibration report
                  if (Config.CalReport) {
                     if (FitSuccess > 0) {
                        CalibrationReport(&Fit, ReportOut, HistName, Settings);
                     } else {
                        ReportOut << endl << "------------------------------------------" << endl << HistName << endl <<
                            "------------------------------------------" << endl << endl;
                        ReportOut << "Fail Fail Fail! The calibration has failed!" << endl;
                     }
                  }
                  // Write .cal file for GRSISpoon
                  if (Config.CalFile) {
                     WriteCalFile(&Fit, CalFileOut, HistName, Settings);
                  }

               } else {
                  if (Config.PrintVerbose) {
                     cout << endl << "Hist " << HistName << " failed to load." << endl;
                  }
               }
            }
         }
      }

   }
   // Now run the fit for the waveform spectrum if required
   // ------------------------------------------------------
   if (Config.CalWave) {
      if (Config.PrintVerbose) {
         cout << "----------------------------------" << endl << "Now fitting Wave Energy Spectra" << endl <<
             "----------------------------------" << endl;
      }

      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            for (Seg = 0; Seg <= SEGS + 1; Seg++) {

               // Set source and histogram name based on seg number
               switch (Seg) {
               case 0:
                  snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02da WaveChg", Clover, Num2Col(Crystal), Seg);
                  HistName = CharBuf;
                  Source = Config.SourceNumCore;
                  break;
               case 9:
                  snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cN%02db WaveChg", Clover, Num2Col(Crystal), 0);
                  HistName = CharBuf;
                  Source = Config.SourceNumCore;
                  break;
               default:
                  snprintf(CharBuf, CHAR_BUFFER_SIZE, "TIG%02d%cP%02dx WaveChg", Clover, Num2Col(Crystal), Seg);
                  HistName = CharBuf;
                  if (Seg < 5) {
                     Source = Config.SourceNumFront;
                  } else {
                     Source = Config.SourceNumBack;
                  }
                  break;
               }

               // Load histogram
               TH1F *Histo = hWaveCharge[Clover - 1][Crystal][Seg];
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
                  Settings.PeakSelect = PeakSelect;
                  Settings.BackupPeakSelect = 0;
                  
                  FitSuccess = FitGammaSpectrum(Histo, &WaveFit, Settings);

                  // Now print reports on results of fits and calibration.
                  if (OUTPUT_GAIN) {
                     if (FitSuccess > 0) {
                        if (FitSuccess < 3 || FORCE_LINEAR) {
                           WaveOut << HistName << ":\t" << WaveFit.LinGainFit[0];
                           WaveOut << "\t" << WaveFit.LinGainFit[1] << endl;
                        } else {
                           WaveOut << HistName << ":\t" << WaveFit.QuadGainFit[0] << "\t";
                           WaveOut << WaveFit.QuadGainFit[1] << "\t" << WaveFit.QuadGainFit[2] << endl;
                        }
                     } else {
                        //WaveOut << HistName << " Fail!!!" << endl;
                     }
                  }
                  if (Config.CalReport) {
                     if (FitSuccess > 0) {
                        WaveReportOut << "test" << endl;
                        CalibrationReport(&WaveFit, WaveReportOut, HistName, Settings);
                     } else {
                        WaveReportOut << endl << "------------------------------------------" << endl << HistName <<
                            endl << "------------------------------------------" << endl << endl;
                        WaveReportOut << "Fail Fail Fail! The calibration has failed!" << endl;
                     }
                  }
               } else {
                  if (Config.PrintVerbose) {
                     cout << endl << "Hist " << HistName << " failed to load." << endl;
                  }
               }
            }
         }
      }

   }

   if (Config.CalEnergy) {
      GainOut.close();
      if (Config.CalReport) {
         ReportOut.close();
      }
      if (Config.CalFile) {
         CalFileOut.close();
      }
   }
   if (Config.CalWave) {
      WaveOut.close();
      if (Config.CalReport) {
         WaveReportOut.close();
      }
   }
   // Write spectra to file again if fits are to be saved
   if (Config.WriteFits) {
      outfile->cd();
      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            for (Seg = 0; Seg <= SEGS + 1; Seg++) {
               dCharge->cd();
               hCharge[Clover - 1][Crystal][Seg]->Write();
               dWaveCharge->cd();
               hWaveCharge[Clover - 1][Crystal][Seg]->Write();
            }
         }
      }
   }



   if (PLOT_CALIB_SUMMARY) {
      cCalib2->cd(1);
      OffsetPlot->Draw();
      cCalib2->Modified();
      cCalib2->Update();
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
   // Write gain histograms to file
   outfile->cd();
   dOther->cd();
   GainPlot->Write();
   OffsetPlot->Write();
   QuadPlot->Write();
   GainHist->Write();
   OffsetHist->Write();
   QuadHist->Write();

   outfile->Close();

}


// 

void ResetTempSpectra()
{
   int Clover, Crystal;
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         hCrystalChargeTemp[Clover - 1][Crystal]->Reset();
      }
   }
}
