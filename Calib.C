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
#include "Utils.h"

// File pointers:
static TFile *outfile = 0;
static TDirectory *dCharge, *dWaveCharge, *dTemp, *dOther = { 0 };

extern TApplication *App;

TCanvas *cWave1;

// Spectra Pointers:
static TH1F *hCharge[CLOVERS][CRYSTALS][SEGS + 2] = { };        // charge from FPGA
static TH1F *hWaveCharge[CLOVERS][CRYSTALS][SEGS + 2] = { };    // charge from waveform

static TH1F *hMidasTime = 0;
static TH1F *hCrystalChargeTemp[CLOVERS][CRYSTALS] = { };       // Only doing these guys for the cores
static TH1F *hCrystalGain[CLOVERS][CRYSTALS] = { };
static TH1F *hCrystalOffset[CLOVERS][CRYSTALS] = { };

static TH1F *WaveHist = 0;

// Functions called from main:   
int InitCalib();
int Calib(std::vector < TTigFragment > &ev);
void FinalCalib();
// Functions called from here:
void ResetTempSpectra();

// Storing run settings
extern RunConfig Config;

int InitCalib()
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
   int Clover, Crystal, Seg;

   if (PLOT_WAVE) {
      cWave1 = new TCanvas();
   }


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
   // Check one source only given
   if (Config.SourceNumCore.size() != 1) {
      cout << "Only one source should be given for calibration of TTree as the sort will sum all input files." << endl;
      return 1;
   }

   if (Config.PrintVerbose) {
      cout << "Searching for core Peaks: " << Config.Sources[Config.
                                                             SourceNumCore[0]][0] << "kev and " <<
          Config.Sources[Config.SourceNumCore[0]]
          [1] << "keV (Ratio " << Config.Sources[Config.SourceNumCore[0]][0] /
          Config.Sources[Config.SourceNumCore[0]][1] << ")" << endl;
      cout << "Searching for front seg peaks: " << Config.Sources[Config.
                                                                  SourceNumFront[0]][0] << "kev and " <<
          Config.Sources[Config.SourceNumFront[0]][1] << "keV (Ratio " << Config.Sources[Config.SourceNumFront[0]][0] /
          Config.Sources[Config.SourceNumFront[0]][1] << ")" << endl;
      cout << "Searching for back seg peaks: " << Config.Sources[Config.
                                                                 SourceNumBack[0]][0] << "kev and " <<
          Config.Sources[Config.SourceNumBack[0]][1] << "keV (Ratio " << Config.Sources[Config.SourceNumBack[0]][0] /
          Config.Sources[Config.SourceNumBack[0]][1] << ")" << endl;
   }


   return 0;
}


int Calib(std::vector < TTigFragment > &ev)
{

   //Variables
   int j, k;
   unsigned int Frag;
   unsigned int Samp, Length;
   int Chan;
   int Crystal, Clover;
   int FitSuccess, CalibSuccess;
   std::string name;

   time_t MidasTime;
   static time_t StartTime;
   static int FirstEvent = 1;
   double RunTimeElapsed;
   static double FitTimeElapsed = 0.0;
   int FileType;
   int FileNum;

   int TimeBin = 0;
   float TB = 0.0;

   float WaveCharge = 0.0;

   if (DEBUG) {
      cout << "--------- New Event ---------" << endl;
   }

   for (Frag = 0; Frag < ev.size(); Frag++) {

      //Get time of first fragment
      if (FirstEvent == 1) {
         StartTime = ev[Frag].MidasTimeStamp;   //ev[Frag].MidasTimeStamp;
         FirstEvent = 0;
         if (Config.PrintBasic) {
            cout << "MIDAS time of first fragment: " << ctime(&StartTime) << endl;
         }
      }
      // Insert test for time earlier than StartTime.....
      if (ev[Frag].MidasTimeStamp < StartTime) {
         StartTime = ev[Frag].MidasTimeStamp;   //ev[i].MidasTimeStamp;
         if (Config.PrintBasic) {
            cout << "Earlier event found! Updating time of first fragment to: " << ctime(&StartTime) << endl;
         }
      }
      //Slave = ((ev[Frag].ChannelAddress & 0x00F00000) >> 20);
      //Port = ((ev[Frag].ChannelAddress & 0x00000F00) >> 8);
      Chan = (ev[Frag].ChannelAddress & 0x000000FF);
      name = ev[Frag].ChannelName;
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
         Length = ev[Frag].wavebuffer.size();
         if (Length > Config.WaveInitialSamples + Config.WaveFinalSamples) {
            //cout << name;// << endl;
            //cout << " samples: " << Length << endl;  
            WaveCharge = CalcWaveCharge(ev[Frag].wavebuffer);
            //cout << "Charge: " << WaveCharge << endl;
            if (PLOT_WAVE) {
               cWave1->cd();
               for (Samp = 0; Samp < Length; Samp++) {
                  WaveHist->SetBinContent(Samp, ev[Frag].wavebuffer.at(Samp));
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
               if (ev[Frag].Charge > 0) {
                  // Increment raw charge spectra
                  if (DEBUG) {
                     cout << "A: Filling " << Clover
                         << ", " << Crystal << ", 0, " << mnemonic.
                         outputsensor << " with charge = " << ev[Frag].Charge << endl;
                  }
                  hCharge[Clover - 1][Crystal][0]->Fill(ev[Frag].Charge);
                  hWaveCharge[Clover - 1][Crystal][0]->Fill(WaveCharge);
                  hCrystalChargeTemp[Clover - 1][Crystal]->Fill(ev[Frag].Charge);
               }
            } else {
               if (Chan == 9) {
                  if (!(mnemonic.outputsensor.compare("b") == 0 || mnemonic.outputsensor.compare("B") == 0)) {  // if this is the primary core output
                     if (Config.PrintBasic) {
                        cout << "Core Channel/Label mismatch" << endl;
                     }
                  }
                  if (ev[Frag].Charge > 0) {
                     // Increment raw charge spectra
                     if (DEBUG) {
                        cout << "B: Filling " << Clover << ", " << Crystal << ", 0, " << mnemonic.outputsensor <<
                            " with charge = " << ev[Frag].Charge << endl;
                     }
                     hCharge[Clover - 1][Crystal][9]->Fill(ev[Frag].Charge);
                     hWaveCharge[Clover - 1][Crystal][9]->Fill(WaveCharge);
                  }
               }
            }
         } else {
            if (mnemonic.segment < 9) {
               if (ev[Frag].Charge > 0) {
                  hCharge[Clover - 1][Crystal][mnemonic.segment]->Fill(ev[Frag].Charge);        // Fill segment spectra
                  hWaveCharge[Clover - 1][Crystal][mnemonic.segment]->Fill(WaveCharge);
               }
            }
         }
         // Now Get time elapsed in run
         MidasTime = ev[Frag].MidasTimeStamp;
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
                  // Configure Fit                  
                  HistoFit Fit;
                  HistoCal Cal;
                  FitSettings Settings = { 0 };
                  FileType = 1; // used to generate histogram name for histo calibration.  Doesn't matter here.  
                  FileNum = 0;  // Used to id source.  Should only be one source type if this function is running.  
                  Config.WriteFits = 0; // Don't want to write fits for these temp spectra.
                  ConfigureEnergyFit(j, k, 0, FileType, FileNum, &Settings);
                  Settings.TempFit = 1;
                  // Perform Fit
                  FitSuccess = FitGammaSpectrum(hCrystalChargeTemp[j - 1][k], &Fit, &Cal, Settings);

                  // Build map of fit results
                  ChannelFitMap ChanFits;
                  for (unsigned int Line = 0; Line < Fit.PeakFits.size(); Line++) {
                     ChanFits.insert(ChannelFitPair(Fit.PeakFits.at(Line).Energy, Fit.PeakFits.at(Line)));
                  }

                  // Clear PeakFits as it will be refilled in CalibrateChannel()
                  Fit.PeakFits.clear();
                  // Calibrate fit map
                  CalibSuccess = CalibrateChannel(ChanFits, Settings, &Fit, &Cal);

                  // Calibration Record
                  if (FitSuccess == 0 && CalibSuccess == 0) {
                     hCrystalGain[j - 1][k]->SetBinContent(TimeBin, Cal.LinGainFit[1]);
                     hCrystalOffset[j - 1][k]->SetBinContent(TimeBin, Cal.LinGainFit[0]);

                  } else {
                     cout << "Calibration of temporary spectra failed!" << endl;
                     continue;
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

   return 0;
}



void FinalCalib()
{

   int Clover = 0;
   int Crystal = 0;
   int Seg = 0;
   string HistName;

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
