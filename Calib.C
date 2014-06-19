// C/C++ libraries:
#include <iostream>
#include <fstream>
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
#include "HistCalib.h"
#include "Utils.h"

// File pointers:
// ---------------
static TFile *outfile = 0;
static TDirectory *dCharge, *dCharge2D, *dWaveCharge, *dWaveCharge2D, *dTemp, *dOther = { 0 };

// ROOT Stuff
//-------------
extern TApplication *App;
TCanvas *cWave1;

// Run settings
extern RunConfig Config;

// Spectra Pointers:
// ----------------
// Charge
static TH1F *hCharge[CLOVERS][CRYSTALS][SEGS + 2] = {0 };        // charge from FPGA
static TH1F *hWaveCharge[CLOVERS][CRYSTALS][SEGS + 2] = {0 };    // charge from waveform
static TH1F *hCrystalChargeTemp[CLOVERS][CRYSTALS] = {0 };       // Only doing these guys for the cores
// 2D Charge
static TH2F *hCoreSegCharge[CLOVERS][CRYSTALS][SEGS] = {0};      // 2D matrix of core charge vs seg charge for low stat seg cal
static TH2F *hCoreSegWaveCharge[CLOVERS][CRYSTALS][SEGS] = {0};
// Calibration 
static TH1F *hCrystalGain[CLOVERS][CRYSTALS] = {0 };    // Histos for recording calibration vs time
static TH1F *hCrystalOffset[CLOVERS][CRYSTALS] = {0 };
// Other
static TH1F *hMidasTime = 0;
static TH1F *hWaveHist = 0;
static TH1F *hWaveChgCrystalFold = 0;
static TH1F *hChgCrystalFold = 0;
static TH2F *hWaveChargeTest = 0;

// Functions
//-------------- 
// called from outside this file:   
int InitCalib();
int Calib(std::vector < TTigFragment > &ev);
void FinalCalib();
// called from here:
void ResetTempSpectra();

int InitCalib()
{
   // Initialise output file                
   std::string tempstring = Config.OutPath + Config.CalOut;
   outfile = new TFile(tempstring.c_str(), "RECREATE");

   dCharge = outfile->mkdir("Charge");
   dCharge2D = outfile->mkdir("Charge2D");
   dWaveCharge = outfile->mkdir("WaveCharge");
   dWaveCharge2D = outfile->mkdir("WaveCharge2D");
   dTemp = outfile->mkdir("Temp");
   dOther = outfile->mkdir("Other");

   char name[CHAR_BUFFER_SIZE], title[CHAR_BUFFER_SIZE];
   int Clover, Crystal, Seg;
   int Scale;

   if (PLOT_WAVE) {
      cWave1 = new TCanvas();
   }

   // Initialise histograms
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         // temp charge histo
         dTemp->cd();
         sprintf(name, "TIG%02d%c Tmp Chg", Clover, Num2Col(Crystal));
         sprintf(title, "TIG%02d%c Temp Core Charge (arb)", Clover, Num2Col(Crystal));
         hCrystalChargeTemp[Clover - 1][Crystal] = new TH1F(name, title, Config.ChargeBins, 0, Config.ChargeMax);
         // histo for record of calibration from temp spectra
         dOther->cd();
         sprintf(name, "TIG%02d%c Gain", Clover, Num2Col(Crystal));
         sprintf(title, "TIG%02d%c Gain vs Time", Clover, Num2Col(Crystal));
         hCrystalGain[Clover - 1][Crystal] = new TH1F(name, title, Config.TimeBins, 0, Config.MaxTime);
         sprintf(name, "TIG%02d%c Offset", Clover, Num2Col(Crystal));
         sprintf(title, "TIG%02d%c Offset vs Time", Clover + 1, Num2Col(Crystal));
         hCrystalOffset[Clover - 1][Crystal] = new TH1F(name, title, Config.TimeBins, 0,Config.MaxTime);
         // Spectra for charge from FPGA
         dCharge->cd();
         Seg = 0;
         sprintf(name, "TIG%02d%cN%02da Chg", Clover, Num2Col(Crystal), Seg);
         sprintf(title, "TIG%02d%cN%02da Charge (arb)", Clover, Num2Col(Crystal), Seg);
         hCharge[Clover - 1][Crystal][Seg] = new TH1F(name, title, Config.ChargeBins, 0, Config.ChargeMax);
         sprintf(name, "TIG%02d%cN%02db Chg", Clover, Num2Col(Crystal), Seg);
         sprintf(title, "TIG%02d%cN%02db Charge (arb)", Clover, Num2Col(Crystal), Seg);
         hCharge[Clover - 1][Crystal][9] = new TH1F(name, title, Config.ChargeBins, 0, Config.ChargeMax);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%cP%02dx Chg", Clover, Num2Col(Crystal), Seg);
            sprintf(title, "TIG%02d%cP%02dx Charge (arb)", Clover, Num2Col(Crystal), Seg);
            hCharge[Clover - 1][Crystal][Seg] = new TH1F(name, title, Config.ChargeBins, 0, Config.ChargeMax);
         }
         // and charge derived from waveform
         dWaveCharge->cd();
         Seg = 0;
         sprintf(name, "TIG%02d%cN%02da WaveChg", Clover, Num2Col(Crystal), Seg);
         sprintf(title, "TIG%02d%cN%02da Waveform Charge (arb)", Clover, Num2Col(Crystal), Seg);
         hWaveCharge[Clover - 1][Crystal][Seg] = new TH1F(name, title, Config.ChargeBins, 0, Config.WaveChargeMax);
         sprintf(name, "TIG%02d%cN%02db WaveChg", Clover, Num2Col(Crystal), Seg);
         sprintf(title, "TIG%02d%cN%02db Waveform Charge (arb)", Clover, Num2Col(Crystal), Seg);
         hWaveCharge[Clover - 1][Crystal][9] = new TH1F(name, title, Config.ChargeBins, 0, Config.WaveChargeMax);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%cP%02dx WaveChg", Clover, Num2Col(Crystal), Seg);
            sprintf(title, "TIG%02d%cP%02dx Waveform Charge (arb)", Clover, Num2Col(Crystal), Seg);
            hWaveCharge[Clover - 1][Crystal][Seg] = new TH1F(name, title, Config.ChargeBins, 0, Config.WaveChargeMax);
         }
         // and 2D charge spectra for low stat seg cal
         if(Config.Cal2D==1) {
            dCharge2D->cd();
            for (Seg = 1; Seg <= SEGS; Seg++) {
               sprintf(name, "TIG%02d%cP%02dx Chg Mat", Clover, Num2Col(Crystal), Seg);
               sprintf(title, "TIG%02d%cP%02dx-TIG%02d%cN00a Charge Matrix  (arb)", Clover, Num2Col(Crystal), Seg, Clover, Num2Col(Crystal));
               hCoreSegCharge[Clover-1][Crystal][Seg-1] = new TH2F(name, title, Config.ChargeBins2D, 0, Config.ChargeMax,  Config.ChargeBins2D, 0, Config.ChargeMax);
            }
            dWaveCharge2D->cd();
            for (Seg = 1; Seg <= SEGS; Seg++) {
               sprintf(name, "TIG%02d%cP%02dx WaveChg Mat", Clover, Num2Col(Crystal), Seg);
               sprintf(title, "TIG%02d%cP%02dx-TIG%02d%cN00a Wave Charge Matrix  (arb)", Clover, Num2Col(Crystal), Seg, Clover, Num2Col(Crystal));
               hCoreSegWaveCharge[Clover-1][Crystal][Seg-1] = new TH2F(name, title, Config.ChargeBins2D, 0, Config.WaveChargeMax,  Config.ChargeBins2D, 0, Config.WaveChargeMax);
            }
         }
      }
   }
   // Time stamp histo
   dOther->cd();
   sprintf(name, "Midas Time");
   sprintf(title, "Midas Timestamps (s)");
   hMidasTime = new TH1F(name, title, Config.TimeBins, 0, Config.MaxTime);

   // Other histograms
   dOther->cd();
   sprintf(name, "Crystal Fold (Chg)");
   sprintf(title, "Crystal Fold (calculated from charge)");
   hChgCrystalFold = new TH1F(name, title, Config.FoldMax, 0, Config.FoldMax);
   sprintf(name, "Crystal Fold (Wave Chg)");
   sprintf(title, "Crystal Fold (calculated from wave charge)");
   hWaveChgCrystalFold = new TH1F(name, title, Config.FoldMax, 0, Config.FoldMax);

   // Wave charge /charge test istogram
      // This spectrum used to check if integration/dispersion are set correct
      // should be a y~=x line, if it is not then something is set wrong and charge thresh will be 
      // different if calculated from waveform or from evaluated energy
   sprintf(name, "Wave Charge Test");
   sprintf(title, "Wave Charge Test");
   Scale = Config.Integration / Config.Dispersion;
   hWaveChargeTest = new TH2F(name, title, Config.ChargeBins2D, 0, Config.ChargeMax/Scale, Config.ChargeBins2D, 0, Config.ChargeMax/Scale);

   if (PLOT_WAVE) {
      sprintf(name, "Wavetemp");
      sprintf(title, "Temporary wave histogram");
      hWaveHist = new TH1F(name, title, Config.WaveformSamples, 0, Config.WaveformSamples);
   }
   // Check one source only given
   if (Config.SourceNumCore.size() != 1) {
      cout << "Exactly one source should be specified for calibration of TTree (sort will sum all input files)." << endl;
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
   int Crystal, Clover, Seg;
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
   
   bool Hits[CLOVERS][CRYSTALS][SEGS+2] = {0};
   bool WaveHits[CLOVERS][CRYSTALS][SEGS+2] = {0};
   int CloverSegFold[CLOVERS] = {0};
   int Charges[CLOVERS][CRYSTALS][SEGS+2] = {0};
   float WaveCharges[CLOVERS][CRYSTALS][SEGS+2] = {0.0};

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

            WaveCharge = CalcWaveCharge(ev[Frag].wavebuffer);

            if (PLOT_WAVE) {
               cWave1->cd();
               for (Samp = 0; Samp < Length; Samp++) {
                  hWaveHist->SetBinContent(Samp, ev[Frag].wavebuffer.at(Samp));
               }

               hWaveHist->Draw();
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
                  if (DEBUG) {
                     cout << "A: Filling " << Clover
                         << ", " << Crystal << ", 0, " << mnemonic.
                         outputsensor << " with charge = " << ev[Frag].Charge << endl;
                  }
                  // Increment histograms
                  hCharge[Clover - 1][Crystal][0]->Fill(ev[Frag].Charge);
                  hWaveCharge[Clover - 1][Crystal][0]->Fill(WaveCharge);
                  hCrystalChargeTemp[Clover - 1][Crystal]->Fill(ev[Frag].Charge);
                  // Store information for use at end of event
                  // Hit records
                  if(TestChargeHit(float(ev[Frag].Charge),Config.Integration,Config.ChargeThresh)) {
                     Hits[Clover - 1][Crystal][0] = 1;
                  }
                  if(TestChargeHit(WaveCharge,1,Config.ChargeThresh)) {
                     WaveHits[Clover - 1][Crystal][0] = 1;
                  }
                  // charge records
                  Charges[Clover - 1][Crystal][0] = ev[Frag].Charge;
                  WaveCharges[Clover - 1][Crystal][0] = WaveCharge;
                  
                  // Increment charge wave charge test
                  if(Clover==1 && Crystal==1) {
                     hWaveChargeTest->Fill(ev[Frag].Charge/Config.Integration,WaveCharge);
                  }
               }
            } else {
               if (Chan == 9) {
                  if (!(mnemonic.outputsensor.compare("b") == 0 || mnemonic.outputsensor.compare("B") == 0)) {  // if this is the primary core output
                     if (Config.PrintBasic) {
                        cout << "Core Channel/Label mismatch" << endl;
                     }
                  }
                  if (ev[Frag].Charge > 0) {
                     if (DEBUG) {
                        cout << "B: Filling " << Clover << ", " << Crystal << ", 0, " << mnemonic.outputsensor <<
                            " with charge = " << ev[Frag].Charge << endl;
                     }
                     // Fill histograms
                     hCharge[Clover - 1][Crystal][9]->Fill(ev[Frag].Charge);
                     hWaveCharge[Clover - 1][Crystal][9]->Fill(WaveCharge);
                     // Store information for use at end of event
                     // Hit records
                     if(TestChargeHit(float(ev[Frag].Charge),Config.Integration,Config.ChargeThresh)) {
                        Hits[Clover - 1][Crystal][9] = 1;
                     }
                     if(TestChargeHit(WaveCharge,1,Config.ChargeThresh)) {
                        WaveHits[Clover - 1][Crystal][9] = 1;
                     }
                     // charge records
                     Charges[Clover - 1][Crystal][9] = ev[Frag].Charge;
                     WaveCharges[Clover - 1][Crystal][9] = WaveCharge;
                     
                  }
               }
            }
         } else {
            if (mnemonic.segment < 9) {
               if (ev[Frag].Charge > 0) {
                  // Fill histograms
                  hCharge[Clover - 1][Crystal][mnemonic.segment]->Fill(ev[Frag].Charge);        // Fill segment spectra
                  hWaveCharge[Clover - 1][Crystal][mnemonic.segment]->Fill(WaveCharge);
                  // Store information for use at end of event
                  // Hit records
                  if(TestChargeHit(float(ev[Frag].Charge),Config.Integration,Config.ChargeThresh)) {
                     Hits[Clover - 1][Crystal][mnemonic.segment] = 1;
                  }
                  if(TestChargeHit(WaveCharge,1,Config.ChargeThresh)) {
                     WaveHits[Clover - 1][Crystal][mnemonic.segment] = 1;
                     // Count segment fold
                     CloverSegFold[Clover-1] += 1;
                  }
                  // charge records
                  Charges[Clover - 1][Crystal][mnemonic.segment] = ev[Frag].Charge;
                  WaveCharges[Clover - 1][Crystal][mnemonic.segment] = WaveCharge;
                  
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
         if ((RunTimeElapsed - FitTimeElapsed >= Config.TimeBinSize) && Config.FitTempSpectra) {

            TB = ((MidasTime - StartTime) / Config.MaxTime);
            TimeBin = TB * Config.TimeBins;

            if (DEBUG) {
               cout << "MT: " << MidasTime << " ST: " << StartTime << " MT: " << Config.MaxTime << " TIME_BINS: " << Config.TimeBins
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
   
   // Now loop hits and increment seg charge - core charge 2D spectra
   if(Config.Cal2D) {
      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            // Check if hit and also if this is the only hit in clover.
            //    (Could do this on a per crystal basis to get more stats but enforcing 1seg hit per clover
            //       means we will certainly not have any problems with crosstalk).
            if(Hits[Clover-1][Crystal][0] == 1 && CloverSegFold[Clover-1] == 1 ) {
               for(Seg=1;Seg<=SEGS;Seg++) {
                  if(Hits[Clover-1][Crystal][Seg] == 1 && WaveHits[Clover-1][Crystal][Seg] == 1) {
                     hCoreSegCharge[Clover-1][Crystal][Seg-1]->Fill(Charges[Clover-1][Crystal][0],Charges[Clover-1][Crystal][Seg]);
                     hCoreSegWaveCharge[Clover-1][Crystal][Seg-1]->Fill(WaveCharges[Clover-1][Crystal][0],WaveCharges[Clover-1][Crystal][Seg]);
                  }
               }
            }     
         }
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
            if(Config.Cal2D && Seg<SEGS) {
               dCharge2D->cd();
               hCoreSegCharge[Clover - 1][Crystal][Seg]->Write();
               dWaveCharge2D->cd();
               hCoreSegWaveCharge[Clover - 1][Crystal][Seg]->Write();
            }
         }
         dTemp->cd();
         hCrystalChargeTemp[Clover - 1][Crystal]->Write();
         dOther->cd();
         hCrystalGain[Clover - 1][Crystal]->Write();
         hCrystalOffset[Clover - 1][Crystal]->Write();
      }
   }
   hMidasTime->Write();
   hWaveChargeTest->Write();
   hChgCrystalFold->Write();
   hWaveChgCrystalFold->Write();
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
