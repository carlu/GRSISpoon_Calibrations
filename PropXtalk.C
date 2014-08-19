// C/C++ libraries:
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
using namespace std;
#include <cstdlib>
#include <math.h>

// ROOT libraries:
#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TCutG.h>
#include <TTreeIndex.h>
#include <TTreePlayer.h>
#include <TChain.h>
#include <TSpectrum.h>
#include <TF1.h>
//#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TApplication.h>

// TriScope libraries
#include "TTigFragment.h"
//#include "TFSPC_Info.h"
//#include "TSharc.h"
//#include "TTigress.h"
//#include "TRf.h"
//#include "TTriFoil.h"

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "PropXtalk.h"
#include "Utils.h"

// stuff
extern TApplication *App;
static TCanvas *cXtalk1;

// File pointers:
static TFile *outfile = 0;
static TDirectory *dRaw, *dWave, *dSum, *dAddBack, *dFold, *dOther = { 0 };

// Spectra pointers  
// for complex spectra naming is h[ITEM][ACTION][RANGE]
//             e.g. SegAddBackCrystal is the add back of segs accross each crystal i.e. total seg energy in that crystal for each event
//                  SegSumCrystal is each individual seg energy which occurs in that crystal
//                  CoreABClover is the total of core energies in that clover for each event

// Raw
static TH1F *hEn[CLOVERS][CRYSTALS][SEGS + 2];  // energy for each individual channel, both cores and segs
static TH1F *hWaveEn[CLOVERS][CRYSTALS][SEGS + 2];      // as above but energy is derived from waveform
static TH1F *hHitPattern;       // record hit counts by TIGRESS DAQ channel numbering
static TH1F *hEHitPattern;      // record above thresh hit counts by TIGRESS DAQ channel numbering
static TH2F *hEnMatrix;         // Channel vs energy matrix for checking calibration
// Sums
static TH1F *hCoreSumTig;       // sum of core energies for array
static TH1F *hCoreSumClover[CLOVERS];   // sum of core energies for each clover
static TH1F *hSegSumClover[CLOVERS];    // sum of seg energies for each clover
static TH1F *hSegSumCrystal[CLOVERS][CRYSTALS]; // sum of seg energies for each crystal
// Addback
static TH1F *hCoreAddBackTig;   // addback of all cores in array for each event
static TH1F *hCoreAddBackClover[CLOVERS];       // addback of all cores in clover for each event
static TH1F *hSegAddBackClover[CLOVERS];        // addback of all segs in clover for each event
//static TH1F *hSegAddBackCrystal[CLOVERS][CRYSTALS];     // addback of all segs in crystal for each event
// Derived
static TH1F *hCloverFoldTig;    // Num clovers hit in array
static TH1F *hCrystalFoldTig;   // Num crystals hit in array
static TH1F *hSegFoldTig;       // Num Segs hit in array
static TH1F *hCrystalFoldClover[CLOVERS];       // Num crystals hit in each clover
static TH1F *hSegFoldClover[CLOVERS];   // Num segs     hit in each clover
static TH1F *hSegFoldCrystal[CLOVERS][CRYSTALS];        // Num segs hit in each Crystal

static TH1F *hFold1CoreEn[CLOVERS][CRYSTALS][SEGS + 1]; // Core energy, for fold one events, for each hit seg
static TH1F *hSegAddBackCloverByFold[CLOVERS][SEGS];    // Sum of segs, over clover, for each event, by seg fold
static TH1F *hSegAddBackCrystalByFold[CLOVERS][CRYSTALS][SEGS]; // Sum od segs, over crystal, for each event, by seg fold

static TH2F *hXTalk[CLOVERS];   // Matrices for dumping XTalk values for inspection
//static TH2F *hXTalkLow[CLOVERS];   // Matrices for dumping XTalk values for inspection

// Storing xtalk
static int XtalkCount[CLOVERS][(SEGS + 2) * CRYSTALS];  // Count corsstalk events for each channel in each clover
static float XTalkFrac[CLOVERS][(SEGS + 2) * CRYSTALS][(SEGS + 2) * CRYSTALS];  // Record crosstalk events for 

// Functions
int InitPropXtalk();
void FinalPropXtalk();
//void SetGains();
float CalibrateEnergy(int Charge, std::vector < float >Coefficients);

void PropXtalk(std::vector < TTigFragment > &ev)
{
   unsigned int Frag;
   int HitClover, HitCrystal, HitSeg;
   int Clover, Crystal, Seg;
   int CloverFoldTig, CrystalFoldTig, SegFoldTig, CrystalFoldClover, SegFoldClover, SegFoldCrystal;
   float CoreABTig, CoreABClover, SegABClover, SegABCrystal;
   float WaveCharge, WaveEnergy;
   char Colour;
   int Chan;
   float En;
   std::string Name;

   float XTalkTemp;
   int XTalkNum, Count;

   unsigned int CalChan;

   int Hits[CLOVERS][CRYSTALS][SEGS + 2] = { {{0}} };
   int CloverCoreFold[CLOVERS] = { 0 }; // + 1 for each core hit in each clover
   int CrystalSegFold[CLOVERS][CRYSTALS] = { {0} };     // +1 for each seg hit in each crystal

   float Energies[CLOVERS][CRYSTALS][SEGS + 2] = { {{0.0}} };
   float WaveEnergies[CLOVERS][CRYSTALS][SEGS + 2] = { {{0.0}} };
   float SegABEn[CLOVERS][CRYSTALS] = { {0.0} };

   int Charges[CLOVERS][CRYSTALS][SEGS + 2] = { {{0}} };

   int CloverHitList[CLOVERS] = { 0 };
   int CloverHitListGood[CLOVERS] = { 0 };

   // Reset stuff
   memset(Hits, 0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(int));
   memset(CloverCoreFold, 0, CLOVERS * sizeof(int));
   memset(CrystalSegFold, 0, CLOVERS * CRYSTALS * sizeof(int));
   memset(Energies, 0.0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(float));
   memset(WaveEnergies, 0.0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(float));
   memset(SegABEn, 0.0, CLOVERS * CRYSTALS * sizeof(float));

   //cout << "Here " << temp++ << " " << ev.size() << endl;

   // --------------------------------------------------------------- //
   // --- First section loops fragments and fills E and hit arrays -- //
   // --------------------------------------------------------------- //

   for (Frag = 0; Frag < ev.size(); Frag++) {
      //Slave = ((ev[Frag].ChannelAddress & 0x00F00000) >> 20);
      //Port = ((ev[Frag].ChannelAddress & 0x00000F00) >> 8);
      Chan = (ev[Frag].ChannelAddress & 0x000000FF);
      Name = ev[Frag].ChannelName;
      WaveCharge = 0.0;
      WaveEnergy = 0.0;
      // Parse name
      Mnemonic mnemonic;
      if (Name.size() >= 10) {
         ParseMnemonic(&Name, &mnemonic);
      } else {
         cout << "Fragment Name Too Short! - This shouldn't happen if the odb is correctly configured!" << endl;
         continue;
      }

      // Fill "Port Hit Pattern"
      hHitPattern->Fill(ev[Frag].ChannelNumber);

      // Get calibrated charge
      if (!Config.HaveAltEnergyCalibration) {
         //cout << "Using standard calibration..." << endl;
         En = ev[Frag].ChargeCal;
      } else {
         int NewCoeffFound = 0;
         std::vector < float >Coefficients;
         for (CalChan = 0; CalChan < Config.EnCalibNames.size(); CalChan++) {
            if (strncmp(Config.EnCalibNames[CalChan].c_str(), Name.c_str(), 9) == 0) { // bug!  this will match the first core
               // name it finds to either a OR b.  Compare 10 chars woud work but then case sensitivity isses on the x/a/b 
               // at the end.  Don't really need second core energy right now so I will come back to this later
               NewCoeffFound = 1;
               break;
            }
         }
         if (NewCoeffFound == 1) {      // If a new set of coeffs was found, then calibrate
            En = CalibrateEnergy(ev[Frag].Charge, Config.EnCalibValues.at(CalChan));
         } else {               // else use the existing calibration
            En = ev[Frag].ChargeCal;
         }
      }

      // Fill "energy hit pattern"
      if (En > Config.EnergyThresh) {
         hEHitPattern->Fill(ev[Frag].ChannelNumber);
      }
      // Get Calibrated "WaveCharge"
      if (ev[Frag].wavebuffer.size() > Config.WaveInitialSamples + Config.WaveFinalSamples) {
         int NewCoeffFound = 0;
         std::vector < float >Coefficients;
         for (CalChan = 0; CalChan < Config.WaveCalibNames.size(); CalChan++) {
            if (strncmp(Config.WaveCalibNames[CalChan].c_str(), Name.c_str(), 9) == 0) {       // bug!  this will match the first core
               // name it finds to either a OR b.  Compare 10 chars woud work but then case sensitivity isses on the x/a/b 
               // at the end.  Don't really need second core energy right now so I will come back to this later
               NewCoeffFound = 1;
               break;
            }
         }
         if (NewCoeffFound == 1) {      // If a new set of coeffs was found, then calibrate
            WaveCharge = CalcWaveCharge(ev[Frag].wavebuffer);
            WaveEnergy = CalibrateWaveEnergy(WaveCharge, Config.WaveCalibValues.at(CalChan));
         }
      } else {
         WaveCharge = 0.0;
         WaveEnergy = 0.0;
      }

      // If TIGRESS
      if (mnemonic.system == "TI") {
         // Determine Crystal
         Colour = mnemonic.arraysubposition.c_str()[0];
         Crystal = Col2Num(Colour);
         if (Crystal == -1) {
            cout << "Bad Colour: " << Colour << endl;
            continue;
         }

         Clover = mnemonic.arrayposition;
         Seg = mnemonic.segment;
         if (Seg == 0) {
            if (Chan == 9) {
               Seg = 9;
            }
         }
         //En  = ev[Frag].ChargeCal;

         //cout << "Here " << temp++ << endl;

         // First record hit pattern, Count Fold, increment raw and sum spectra
         Energies[Clover - 1][Crystal][Seg] = En;
         WaveEnergies[Clover - 1][Crystal][Seg] = WaveEnergy;
         Charges[Clover - 1][Crystal][Seg] = ev[Frag].Charge;
         CloverHitList[Clover - 1] |= 1;
         if (En > Config.EnergyThresh) {
            CloverHitListGood[Clover - 1] |= 1;
         }


         switch (Seg) {
         case 0:
            if (En > Config.EnergyThresh) {
               Hits[Clover - 1][Crystal][Seg] += 1;     // This should only ever be 1, if 2 or more then we have an event building problem
               CloverCoreFold[Clover - 1] += 1;
               hCoreSumTig->Fill(En);
               hCoreSumClover[Clover - 1]->Fill(En);
               hEn[Clover - 1][Crystal][Seg]->Fill(En);
               hEnMatrix->Fill(ev[Frag].ChannelNumber, En);
            }
            if (WaveEnergy > Config.EnergyThresh) {
               hWaveEn[Clover - 1][Crystal][Seg]->Fill(WaveEnergy);
            }
            break;
         case 9:
            if (En > Config.EnergyThresh) {
               Hits[Clover - 1][Crystal][Seg] += 1;
               hEn[Clover - 1][Crystal][Seg]->Fill(En);
               hEnMatrix->Fill(ev[Frag].ChannelNumber, En);
            }
            if (WaveEnergy > Config.EnergyThresh) {
               hWaveEn[Clover - 1][Crystal][Seg]->Fill(WaveEnergy);
            }
            break;
         default:              // Should catch segs only
            if (Seg < 1 || Seg > 8) {
               cout << "Seg numbering problem!!" << endl;
            }
            if (En > Config.EnergyThresh) {
               Hits[Clover - 1][Crystal][Seg] += 1;
               SegABEn[Clover - 1][Crystal] += En;
               CrystalSegFold[Clover - 1][Crystal] += 1;
               hEn[Clover - 1][Crystal][Seg]->Fill(En);
               hEnMatrix->Fill(ev[Frag].ChannelNumber, En);
               hSegSumClover[Clover - 1]->Fill(En);
               hSegSumCrystal[Clover - 1][Crystal]->Fill(En);
            }
            if (WaveEnergy > Config.EnergyThresh) {
               hWaveEn[Clover - 1][Crystal][Seg]->Fill(WaveEnergy);
            }
            break;

         }
      }
   }


   // --------------------------------------------------------------- //
   // --- Debug section.  Checks hits, energies, add-back          -- //
   // --------------------------------------------------------------- //

   if (DEBUG_HITS) {
      cout << endl << "-------------------" << endl << " Clover Hits" << endl << "------------------" << endl;
      cout << "Hit:\t";
      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         if (CloverHitList[Clover - 1] == 1) {
            cout << "1\t";
         } else {
            cout << "0\t";
         }
      }
      cout << endl;
      cout << "GoodE:\t";
      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         if (CloverHitListGood[Clover - 1] == 1) {
            cout << "1\t";
         } else {
            cout << "0\t";
         }
      }
      cout << endl;

      cout << endl << "-------------------" << endl << "Segment Hits and Energies" << endl << "------------------" <<
          endl;
      for (Clover = 1; Clover <= CLOVERS; Clover++) {
         if (CloverHitList[Clover - 1] == 1) {
            cout << endl;
            for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
               cout << "Clover " << Clover << " Crystal " << Crystal << endl;
               for (Seg = 0; Seg < SEGS + 2; Seg++) {
                  cout << Charges[Clover - 1][Crystal][Seg] << "\t\t";
               }
               cout << endl;
               for (Seg = 0; Seg < SEGS + 2; Seg++) {
                  cout << Energies[Clover - 1][Crystal][Seg] << "\t";
               }
               cout << endl;
               for (Seg = 0; Seg < SEGS + 2; Seg++) {
                  cout << Hits[Clover - 1][Crystal][Seg] << "\t\t";
               }
               cout << endl;
            }
         }
      }
      cout << endl;
   }                            // End of hit debugging section




   //cout << "Hello" << endl;
   // Loop hit patterns, calculate folds, fill addback and fold spectra
   CrystalFoldTig = 0;
   SegFoldTig = 0;
   CloverFoldTig = 0;
   CoreABTig = 0.0;
   HitSeg = 0;
   HitClover = 0; 
   HitCrystal = 0;

   // --------------------------------------------------------------- //
   // --- Final section loops detector elements, calculates        -- //
   // ---    derrived quantities like fold, applies gates,         -- //
   // ---    calculates crosstalk                                  -- //
   // --------------------------------------------------------------- //

   for (Clover = 1; Clover <= CLOVERS; Clover++) {

      //cout << "Clover " << Clover << endl;

      SegFoldClover = 0;
      CoreABClover = 0.0;
      SegABClover = 0.0;

      if (CloverCoreFold[Clover - 1] > 0) {
         CloverFoldTig += 1;
      }
      CrystalFoldClover = 0;
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         SegFoldCrystal = 0;
         SegABCrystal = 0.0;
         if (Hits[Clover - 1][Crystal][0] > 0) {
            CrystalFoldClover += 1;
            CrystalFoldTig += 1;
            CoreABTig += Energies[Clover - 1][Crystal][0];
            CoreABClover += Energies[Clover - 1][Crystal][0];
            HitClover = Clover;
            HitCrystal = Crystal;
         }
         for (Seg = 1; Seg <= SEGS; Seg++) {
            if (Hits[Clover - 1][Crystal][Seg] > 0) {
               SegFoldTig += 1;
               SegFoldClover += 1;
               SegFoldCrystal += 1;
               SegABClover += Energies[Clover - 1][Crystal][Seg];
               SegABCrystal += Energies[Clover - 1][Crystal][Seg];
               HitSeg = Seg;
            }
         }
         if (SegFoldCrystal > 0 || Hits[Clover - 1][Crystal][0] > 0) {      // If hit here in this crystal core OR any seg, then record seg fold.
            hSegFoldCrystal[Clover - 1][Crystal]->Fill(SegFoldCrystal);
            hSegAddBackCrystalByFold[Clover - 1][Crystal][0]->Fill(SegABEn[Clover - 1][Crystal]);
            hSegAddBackCrystalByFold[Clover - 1][Crystal][SegFoldCrystal]->Fill(SegABEn[Clover - 1][Crystal]);
         }
         //hSegAddBackCrystal[Clover - 1][Crystal]->Fill(SegABCrystal);

      }

      //cout << "CrystalFoldClover: " << CrystalFoldClover << " SegFoldClover: " << SegFoldClover << endl;
      //cout << "CoreABClover: " << CoreABClover << " SegABClover: " << SegABClover << endl;

      if (CrystalFoldClover > 0 || SegFoldClover > 0) { // If hit in core OR seg of this clover, then increment fold
         hCrystalFoldClover[Clover - 1]->Fill(CrystalFoldClover);
         hSegFoldClover[Clover - 1]->Fill(SegFoldClover);
         hCoreAddBackClover[Clover - 1]->Fill(CoreABClover);
         hSegAddBackClover[Clover - 1]->Fill(SegABClover);
      }
      if (SegFoldClover == 1) {
         // Xtalk calculation here
         // checks - One seg hit  
         //        - CoreE = SegE
         //        - CoreABCloverE = CoreE 
         // As this is fold1, should be able to use HitClover, HitCrystal and HitSeg.

         // Energy Gate
         if (Energies[Clover - 1][HitCrystal][HitSeg] > 100.0) {
            // Check CoreE = SegE, CoreABCloverE = CoreE    
            //cout << "Cl: " << Clover << " Cr: " << HitCrystal << " HS: " << HitSeg << " En: " << Energies[Clover - 1][HitCrystal][HitSeg] << " CoreEn: " << Energies[Clover - 1][HitCrystal][0] << endl;
            //cout << "Cl: " << Clover << " Cr: " << HitCrystal <<" HS: " << HitSeg << " WEn: " << WaveEnergies[Clover - 1][HitCrystal][HitSeg] << " CorewEn: " << WaveEnergies[Clover - 1][HitCrystal][0] << endl << endl;
            // Count Events

            Count = XtalkCount[Clover - 1][(HitCrystal * 4) + HitSeg];

            // Loop and record crosstalk
            for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
               for (Seg = 0; Seg < SEGS + 2; Seg++) {
                  XTalkTemp = WaveEnergies[Clover - 1][Crystal][Seg] / Energies[Clover - 1][HitCrystal][HitSeg];
                  XTalkNum = (HitCrystal * (SEGS + 2)) + HitSeg;
                  XTalkFrac[Clover - 1][XTalkNum][(Crystal * (SEGS + 2)) + Seg] =
                      ((XTalkFrac[Clover - 1][XTalkNum][(Crystal * (SEGS + 2)) + Seg] * Count) + XTalkTemp) / (Count + 1);
               }
            }
            XtalkCount[Clover - 1][(HitCrystal * 4) + HitSeg] += 1;

         }

      }
   }

   hCloverFoldTig->Fill(CloverFoldTig);
   hCrystalFoldTig->Fill(CrystalFoldTig);
   hSegFoldTig->Fill(SegFoldTig);
   hCoreAddBackTig->Fill(CoreABTig);

}


int InitPropXtalk()
{
   char Colours[] = "BGRW";
   char name[512], title[512];
   int Clover, Crystal, Seg, Fold;

   // Initialise ROOT output file   
   std::string tempstring = Config.OutPath + Config.PropOut;
   outfile = new TFile(tempstring.c_str(), "RECREATE");

   dRaw = outfile->mkdir("Raw");
   dSum = outfile->mkdir("Sum");
   dAddBack = outfile->mkdir("Add-Back");
   dFold = outfile->mkdir("Fold");
   dOther = outfile->mkdir("Other");
   dWave = outfile->mkdir("Wave");
   
   // Initialise spectra
   // Raw
   dRaw->cd();
   // Create Spectra
   sprintf(name, "Hit Pattern");
   sprintf(title, "TIGRESS Port Hit Pattern");
   hHitPattern = new TH1F(name, title, 2000, 0, 2000);
   sprintf(name, "EHit Pattern");
   sprintf(title, "TIGRESS Energy Hit Pattern");
   hEHitPattern = new TH1F(name, title, 2000, 0, 2000);
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         Seg = 0;
         sprintf(name, "TIG%02d%cN%02dA Core En", Clover , Colours[Crystal], Seg);
         sprintf(title, "TIG%02d%cN%02dA Core A Energy (keV)", Clover, Colours[Crystal], Seg);
         hEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%cP%02dx Seg En", Clover , Colours[Crystal], Seg);
            sprintf(title, "TIG%02d%cP%02dx Seg Energy (keV)", Clover, Colours[Crystal], Seg);
            hEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         }
         Seg = SEGS + 1;
         sprintf(name, "TIG%02d%cN%02dB Core En", Clover, Colours[Crystal], 0);
         sprintf(title, "TIG%02d%cN%02dB Core B Energy (keV)", Clover, Colours[Crystal], 0);
         hEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }
   sprintf(name, "En v Channel");
   sprintf(title, "TIGRESS DAQ Channel vs Calibrated Energy (keV)");
   hEnMatrix = new TH2F(name, title, 1000, 0, 1000, 2000, 0, 2000);
   // Waveform energy
   dWave->cd();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         Seg = 0;
         sprintf(name, "TIG%02d%cN%02dA Core WaveEn", Clover, Colours[Crystal], Seg);
         sprintf(title, "TIG%02d%cN%02dA Core A Waveform Energy (keV)", Clover, Colours[Crystal], Seg);
         hWaveEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%cP%02dx Seg WaveEn", Clover, Colours[Crystal], Seg);
            sprintf(title, "TIG%02d%cP%02dx Seg Waveform Energy (keV)", Clover, Colours[Crystal], Seg);
            hWaveEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         }
         Seg = SEGS + 1;
         sprintf(name, "TIG%02d%cN%02dB Core WaveEn", Clover, Colours[Crystal], 0);
         sprintf(title, "TIG%02d%cN%02dB Core B Waveform Energy (keV)", Clover, Colours[Crystal], 0);
         hWaveEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }

   //Sum
   dSum->cd();
   sprintf(name, "TIGRESS Core Sum");
   sprintf(title, "TIGRESS Core Sum Energy (keV)");
   hCoreSumTig = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      sprintf(name, "TIG%02d Core Sum", Clover);
      sprintf(title, "TIG%02d Core Sum Energy (keV)", Clover);
      hCoreSumClover[Clover - 1] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      sprintf(name, "TIG%02d Seg Sum", Clover);
      sprintf(title, "TIG%02d Segment Sum Energy (keV)", Clover);
      hSegSumClover[Clover - 1] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         sprintf(name, "TIG%02d%c Seg Sum", Clover, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Segment Sum Energy (keV)", Clover, Colours[Crystal]);
         hSegSumCrystal[Clover - 1][Crystal] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }
   //Addback
   dAddBack->cd();
   sprintf(name, "TIGRESS Core AB");
   sprintf(title, "TIGRESS All Core Add-Back Energy (keV)");
   hCoreAddBackTig = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      sprintf(name, "TIG%02d Core AB", Clover);
      sprintf(title, "TIG%02d Core Add-Back Energy (keV)", Clover);
      hCoreAddBackClover[Clover - 1] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      sprintf(name, "TIG%02d Seg AB", Clover);
      sprintf(title, "TIG%02d Segment Add-Back Energy (keV)", Clover);
      hSegAddBackClover[Clover - 1] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      //Seg AB Crystal???
   }
   //Derived
   dFold->cd();
   sprintf(name, "TIG Clover Fold");
   sprintf(title, "TIGRESS Clover Fold");
   hCloverFoldTig = new TH1F(name, title, Config.FoldMax, 0, Config.FoldMax);
   sprintf(name, "TIG Crys Fold");
   sprintf(title, "TIGRESS Crystal Fold");
   hCrystalFoldTig = new TH1F(name, title, Config.FoldMax, 0, Config.FoldMax);
   sprintf(name, "TIG Seg Fold");
   sprintf(title, "TIGRESS Segment Fold");
   hSegFoldTig = new TH1F(name, title, Config.FoldMax, 0, Config.FoldMax);
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      sprintf(name, "TIG%02d Crys Fold", Clover);
      sprintf(title, "TIG%02d Crystal Fold", Clover);
      hCrystalFoldClover[Clover - 1] = new TH1F(name, title, Config.FoldMax, 0, Config.FoldMax);
      sprintf(name, "TIG%02d Seg Fold", Clover);
      sprintf(title, "TIG%02d Segment Fold", Clover);
      hSegFoldClover[Clover - 1] = new TH1F(name, title, Config.FoldMax, 0, Config.FoldMax);
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         sprintf(name, "TIG%02d%c Seg Fold", Clover, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Segment Fold", Clover, Colours[Crystal]);
         hSegFoldCrystal[Clover - 1][Crystal] = new TH1F(name, title, Config.FoldMax, 0, Config.FoldMax);
      }
   }
   dOther->cd();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         sprintf(name, "TIG%02d%c Core En F1", Clover, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Core Energy SegFold1 Any Seg (keV)", Clover, Colours[Crystal]);
         hFold1CoreEn[Clover - 1][Crystal][0] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%c Core En F1 Seg%02d", Clover, Colours[Crystal], Seg);
            sprintf(title, "TIG%02d%c Core Energy SegFold1 Segment%02d (keV)", Clover, Colours[Crystal], Seg);
            hFold1CoreEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         }
      }
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         Fold = 0;
         sprintf(name, "TIG%02d%c Seg AB", Clover, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Segment Add-Back Energy Any Fold (keV)", Clover, Colours[Crystal]);
         hSegAddBackCrystalByFold[Clover - 1][Crystal][Fold] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         for (Fold = 1; Fold < SEGS; Fold++) {
            sprintf(name, "TIG%02d%c Seg AB F%01d", Clover, Colours[Crystal], Fold);
            sprintf(title, "TIG%02d%c Segment Add-Back Energy Fold %01d (keV)", Clover, Colours[Crystal], Fold);
            hSegAddBackCrystalByFold[Clover - 1][Crystal][Fold] =
                new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         }
      }

      sprintf(name, "TIG%02d Seg AB F0", Clover);
      sprintf(title, "TIG%02d Segment Add-Back Energy Any Fold (keV)", Clover);
      hSegAddBackCloverByFold[Clover - 1][0] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      for (Fold = 1; Fold < SEGS; Fold++) {
         sprintf(name, "TIG%02d Seg AB F%01d", Clover, Fold);
         sprintf(title, "TIG%02d Segment Add-Back Energy Fold %01d (keV)", Clover, Fold);
         hSegAddBackCloverByFold[Clover - 1][Fold] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }

   // 2D crosstalk matrices 
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      sprintf(name, "TIG%02d Xtalk", Clover);
      sprintf(title, "TIG%02d Fractional Crosstalk", Clover);
      hXTalk[Clover - 1] =
          new TH2F(name, title, (CRYSTALS * (SEGS + 2)), 0, (CRYSTALS * (SEGS + 2)), (CRYSTALS * (SEGS + 2)), 0,
                   (CRYSTALS * (SEGS + 2)));
      printf(name, "TIG%02d XtalkLow", Clover);
      sprintf(title, "TIG%02d Fractional Crosstalk (Low values)", Clover);
      //hXTalkLow[Clover - 1] =
      //  new TH2F(name, title, (CRYSTALS * (SEGS + 2)), 0, (CRYSTALS * (SEGS + 2)), (CRYSTALS * (SEGS + 2)), 0,
      //         (CRYSTALS * (SEGS + 2)));
   }
   return 0;
}


void FinalPropXtalk()
{

   outfile->cd();
   if (PLOT_ON == 1) {
      cXtalk1 = new TCanvas();
      cXtalk1->cd();
   }
   int Clover, Crystal, Seg, Fold;
   int HitSegment, OtherSegment;

   // Write crosstalk values out to text file and fill 2d hist
   ofstream XTalkOut;
   std::string tempstring = Config.OutPath + Config.PropTxtOut;
   XTalkOut.open(tempstring.c_str());
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      XTalkOut << endl << "------------------" << endl << " Clover " << Clover
           << endl << "------------------" << endl << endl;
      for (HitSegment = 0; HitSegment < ((SEGS + 2) * CRYSTALS); HitSegment++) {
         for (OtherSegment = 0; OtherSegment < ((SEGS + 2) * CRYSTALS); OtherSegment++) {
            XTalkOut << XTalkFrac[Clover - 1][HitSegment][OtherSegment] << " ";
            hXTalk[Clover - 1]->SetBinContent(HitSegment + 1, OtherSegment + 1,
                                          XTalkFrac[Clover - 1][HitSegment][OtherSegment]);
            //hXTalkLow[Clover - 1]->SetBinContent(HitSegment+1, OtherSegment+1, XTalkFrac[Clover - 1][HitSegment][OtherSegment]);
         }
         XTalkOut << endl;
      }
   }
   XTalkOut.close();

   // Write spectra to file
   // Raw
   outfile->cd();
   dRaw->cd();
   hHitPattern->Write();
   hEHitPattern->Write();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg < SEGS + 2; Seg++) {
            hEn[Clover - 1][Crystal][Seg]->Write();
         }
      }
   }
   hEnMatrix->Write();
   // Waveform energy
   dWave->cd();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg < SEGS + 2; Seg++) {
            hWaveEn[Clover - 1][Crystal][Seg]->Write();
         }
      }
   }
   // Sums
   dSum->cd();
   hCoreSumTig->Write();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      hCoreSumClover[Clover - 1]->Write();
      hSegSumClover[Clover - 1]->Write();
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         hSegSumCrystal[Clover - 1][Crystal]->Write();
      }
   }

   // Add-Back
   dAddBack->cd();
   hCoreAddBackTig->Write();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      hCoreAddBackClover[Clover - 1]->Write();
      hSegAddBackClover[Clover - 1]->Write();
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         //hSegAddBackCrystal[Clover - 1][Crystal]->Write();
      }
   }

   // Derived
   dFold->cd();
   hCloverFoldTig->Write();
   hCrystalFoldTig->Write();
   hSegFoldTig->Write();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      hCrystalFoldClover[Clover - 1]->Write();      // Num crystals hit in each clover
      hSegFoldClover[Clover - 1]->Write();
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         hSegFoldCrystal[Clover - 1][Crystal]->Write();
      }
   }
   dOther->cd();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg < SEGS + 1; Seg++) {
            hFold1CoreEn[Clover - 1][Crystal][Seg]->Write();
         }
      }
   }

   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Fold = 0; Fold < SEGS; Fold++) {
         hSegAddBackCloverByFold[Clover - 1][Fold]->Write();
      }
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Fold = 0; Fold < SEGS; Fold++) {
            hSegAddBackCrystalByFold[Clover - 1][Crystal][Fold]->Write();
         }
      }
   }

   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      hXTalk[Clover - 1]->Write();
      //hXTalkLow[Clover - 1]->Write();
   }

   outfile->Close();

}
