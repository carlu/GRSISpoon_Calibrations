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
#include "Main.h"
#include "PropXtalk.h"
//#include "Gains.h"

// Gain arrays
//#include "Gains.h"

// stuff
extern TApplication *App;
static TCanvas *cXtalk1;

// File pointers:
static TFile *outfile = 0;
static TDirectory *dRaw, *dSum, *dAddBack, *dFold, *dOther = { 0 };

// Spectra pointers  
// for complex spectra naming is h[ITEM][ACTION][RANGE]
//             e.g. SegAddBackCrystal is the add back of segs accross each crystal i.e. total seg energy in that crystal for each event
//                  SegSumCrystal is each individual seg energy which occurs in that crystal
//                  CoreABClover is the total of core energies in that clover for each event

// Raw
static TH1F *hEn[CLOVERS][CRYSTALS][SEGS + 2];  // energy for each individual channel, both cores and segs
static TH1F *hHitPattern;       // record hit counts by TIGRESS DAQ channel numbering
static TH1F *hEHitPattern;      // record above thresh hit counts by TIGRESS DAQ channel numbering
// Sums
static TH1F *hCoreSumTig;       // sum of core energies for array
static TH1F *hCoreSumClover[CLOVERS];   // sum of core energies for each clover
static TH1F *hSegSumClover[CLOVERS];    // sum of seg energies for each clover
static TH1F *hSegSumCrystal[CLOVERS][CRYSTALS]; // sum of seg energies for each crystal
// Addback
static TH1F *hCoreAddBackTig;   // addback of all cores in array for each event
static TH1F *hCoreAddBackClover[CLOVERS];       // addback of all cores in clover for each event
static TH1F *hSegAddBackClover[CLOVERS];        // addback of all segs in clover for each event
static TH1F *hSegAddBackCrystal[CLOVERS][CRYSTALS];     // addback of all segs in crystal for each event
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

// Storing xtalk
static int XtalkCount[CLOVERS][(SEGS + 2) * CRYSTALS];  // Count corsstalk events for each channel in each clover
static float XTalkFrac[CLOVERS][(SEGS + 2) * CRYSTALS][(SEGS + 2) * CRYSTALS];  // Record crosstalk events for 

// Functions
void InitPropXtalk();
void FinalPropXtalk();
//void SetGains();
float CalibrateEnergy(int Charge, std::vector < float >Coefficients);

void PropXtalk(std::vector < TTigFragment > &ev)
{

   int temp = 0;
   //cout << "Here " << temp++ << endl;

   int i;
   int HitClover, HitCrystal, HitSeg;
   int Clover, Crystal, Seg, Fold;
   int CloverFoldTig, CrystalFoldTig, SegFoldTig, CrystalFoldClover, SegFoldClover, SegFoldCrystal;
   float CoreABTig, CoreABClover, SegABClover, SegABCrystal;
   char Colour;
   int Slave, Port, Chan;
   float En;
   std::string Name;

   float XTalkTemp;
   int XTalkNum, Count;

   int CalChan;

   int Hits[CLOVERS][CRYSTALS][SEGS + 2] = { 0 };
   int CloverCoreFold[CLOVERS] = { 0 }; // + 1 for each core hit in each clover
   int CrystalSegFold[CLOVERS][CRYSTALS] = { 0 };       // +1 for each seg hit in each crystal

   float Energies[CLOVERS][CRYSTALS][SEGS + 2] = { 0.0 };
   float SegABEn[CLOVERS][CRYSTALS] = { 0.0 };

   int Charges[CLOVERS][CRYSTALS][SEGS + 2] = { 0 };

   int CloverHitList[CLOVERS] = { 0 };
   int CloverHitListGood[CLOVERS] = { 0 };

   // Reset stuff
   memset(Hits, 0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(int));
   memset(CloverCoreFold, 0, CLOVERS * sizeof(int));
   memset(CrystalSegFold, 0, CLOVERS * CRYSTALS * sizeof(int));
   memset(Energies, 0.0, CLOVERS * CRYSTALS * (SEGS + 2) * sizeof(float));
   memset(SegABEn, 0.0, CLOVERS * CRYSTALS * sizeof(float));

   //cout << "Here " << temp++ << " " << ev.size() << endl;

   // --------------------------------------------------------------- //
   // --- First section loops fragments and fills E and hit arrays -- //
   // --------------------------------------------------------------- //

   for (i = 0; i < ev.size(); i++) {
      Slave = ((ev[i].ChannelAddress & 0x00F00000) >> 20);
      Port = ((ev[i].ChannelAddress & 0x00000F00) >> 8);
      Chan = (ev[i].ChannelAddress & 0x000000FF);
      Name = ev[i].ChannelName;

      // Parse name
      Mnemonic mnemonic;
      if (Name.size() >= 10) {
         ParseMnemonic(&Name, &mnemonic);
      } else {
         cout << "Fragment Name Too Short! - This shouldn't happen if the odb is correctly configured!" << endl;
         continue;
      }

      // Fill "Port Hit Pattern"
      hHitPattern->Fill(ev[i].ChannelNumber);


      // Get calibrated charge

      if (!USE_ALT_CALIB) {
         //cout << "Using standard calibration..." << endl;
         En = ev[i].ChargeCal;
      } else {
         int NewCoeffFound = 0;
         std::vector < float >Coefficients;
         for (CalChan = 0; CalChan < CalibNames.size(); CalChan++) {
            if (strncmp(CalibNames[CalChan].c_str(), Name.c_str(), 9) == 0) {
               NewCoeffFound = 1;
               break;
            }
         }
         if (NewCoeffFound == 1) {      // If a new set of coeffs was found, then calibrate
            //En = CalibrateEnergy(ev[i].Charge,Coefficients);
            En = CalibrateEnergy(ev[i].Charge, CalibValues.at(CalChan));
         } else {               // else use the existing calibration
            En = ev[i].ChargeCal;
         }
      }


      // Fill "energy hit pattern"
      if (En > EN_THRESH) {
         hEHitPattern->Fill(ev[i].ChannelNumber);
      }
      //cout << "Here " << temp++ << endl;

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
         //En  = ev[i].ChargeCal;

         //cout << "Here " << temp++ << endl;

         // First record hit pattern, Count Fold, increment raw and sum spectra
         Energies[Clover - 1][Crystal][Seg] = En;
         Charges[Clover - 1][Crystal][Seg] = ev[i].Charge;
         CloverHitList[Clover - 1] |= 1;
         if (En > EN_THRESH) {
            CloverHitListGood[Clover - 1] |= 1;
         }

         switch (Seg) {
         case 0:
            if (En > EN_THRESH) {
               Hits[Clover - 1][Crystal][Seg] += 1;     // This should only ever be 1, if 2 or more then we have an event building problem
               CloverCoreFold[Clover - 1] += 1;
               hCoreSumTig->Fill(En);
               hCoreSumClover[Clover - 1]->Fill(En);
               hEn[Clover - 1][Crystal][Seg]->Fill(En);
            }
            break;
         case 9:
            if (En > EN_THRESH) {
               Hits[Clover - 1][Crystal][Seg] += 1;
               hEn[Clover - 1][Crystal][Seg]->Fill(En);
            }
            break;
         default:              // Should catch segs only
            if (Seg < 1 || Seg > 8) {
               cout << "Seg numbering problem!!" << endl;
            }
            if (En > EN_THRESH) {
               Hits[Clover - 1][Crystal][Seg] += 1;
               SegABEn[Clover - 1][Crystal] += En;
               CrystalSegFold[Clover - 1][Crystal] += 1;
               hEn[Clover - 1][Crystal][Seg]->Fill(En);
               hSegSumClover[Clover - 1]->Fill(En);
               hSegSumCrystal[Clover - 1][Crystal]->Fill(En);
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
      for (Clover = 0; Clover < CLOVERS; Clover++) {
         if (CloverHitList[Clover] == 1) {
            cout << "1\t";
         } else {
            cout << "0\t";
         }
      }
      cout << endl;
      cout << "GoodE:\t";
      for (Clover = 0; Clover < CLOVERS; Clover++) {
         if (CloverHitListGood[Clover] == 1) {
            cout << "1\t";
         } else {
            cout << "0\t";
         }
      }
      cout << endl;

      cout << endl << "-------------------" << endl << "Segment Hits and Energies" << endl << "------------------" <<
          endl;
      for (Clover = 0; Clover < CLOVERS; Clover++) {
         if (CloverHitList[Clover] == 1) {
            cout << endl;
            for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
               cout << "Clover " << Clover + 1 << " Crystal " << Crystal << endl;
               for (Seg = 0; Seg < SEGS + 2; Seg++) {
                  cout << Charges[Clover][Crystal][Seg] << "\t\t";
               }
               cout << endl;
               for (Seg = 0; Seg < SEGS + 2; Seg++) {
                  cout << Energies[Clover][Crystal][Seg] << "\t";
               }
               cout << endl;
               for (Seg = 0; Seg < SEGS + 2; Seg++) {
                  cout << Hits[Clover][Crystal][Seg] << "\t\t";
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


   // --------------------------------------------------------------- //
   // --- Final section loops detector elements, calculates        -- //
   // ---    derrived quantities like fold, applies gates,         -- //
   // ---    calculates crosstalk                                  -- //
   // --------------------------------------------------------------- //

   for (Clover = 0; Clover < CLOVERS; Clover++) {

      //cout << "Clover " << Clover << endl;

      SegFoldClover = 0;
      CoreABClover = 0.0;
      SegABClover = 0.0;

      if (CloverCoreFold[Clover] > 0) {
         CloverFoldTig += 1;
      }
      CrystalFoldClover = 0;
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         SegFoldCrystal = 0;
         SegABCrystal = 0.0;
         if (Hits[Clover][Crystal][0] > 0) {
            CrystalFoldClover += 1;
            CrystalFoldTig += 1;
            CoreABTig += Energies[Clover][Crystal][0];
            CoreABClover += Energies[Clover][Crystal][0];
            HitClover = Clover;
            HitCrystal = Crystal;
         }
         for (Seg = 1; Seg <= SEGS + 1; Seg++) {
            if (Hits[Clover][Crystal][Seg] > 0) {
               SegFoldTig += 1;
               SegFoldClover += 1;
               SegFoldCrystal += 1;
               SegABClover += Energies[Clover][Crystal][Seg];
               SegABCrystal += Energies[Clover][Crystal][Seg];
               HitSeg = Seg;
            }
         }
         if (SegFoldCrystal > 0 || Hits[Clover][Crystal][0] > 0) {      // If hit here in this crystal core OR any seg, then record seg fold.
            hSegFoldCrystal[Clover][Crystal]->Fill(SegFoldCrystal);
            hSegAddBackCrystalByFold[Clover][Crystal][0]->Fill(SegABEn[Clover][Crystal]);
            hSegAddBackCrystalByFold[Clover][Crystal][SegFoldCrystal]->Fill(SegABEn[Clover][Crystal]);
         }
         //hSegAddBackCrystal[Clover][Crystal]->Fill(SegABCrystal);

      }

      //cout << "CrystalFoldClover: " << CrystalFoldClover << " SegFoldClover: " << SegFoldClover << endl;
      //cout << "CoreABClover: " << CoreABClover << " SegABClover: " << SegABClover << endl;

      if (CrystalFoldClover > 0 || SegFoldClover > 0) { // If hit in core OR seg of this clover, then increment fold
         hCrystalFoldClover[Clover]->Fill(CrystalFoldClover);
         hSegFoldClover[Clover]->Fill(SegFoldClover);
         hCoreAddBackClover[Clover]->Fill(CoreABClover);
         hSegAddBackClover[Clover]->Fill(SegABClover);
      }
      if (SegFoldClover == 1) {
         // Xtalk calculation here
         // checks - One seg hit  
         //        - CoreE = SegE
         //        - CoreABCloverE = CoreE    
         // As this is fold1, should be able to use HitClover, HitCrystal and HitSeg.

         // Energy Gate
         if (Energies[Clover][HitCrystal][HitSeg] > 100.0) {
            // Check CoreE = SegE, CoreABCloverE = CoreE    

            // Count Events

            Count = XtalkCount[Clover][(HitCrystal * 4) + HitSeg];

            // Loop and record crosstalk
            for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
               for (Seg = 0; Seg < SEGS + 2; Seg++) {
                  XTalkTemp = Energies[Clover][Crystal][Seg] / Energies[Clover][HitCrystal][HitSeg];
                  XTalkNum = (HitCrystal * (SEGS + 2)) + HitSeg;
                  XTalkFrac[Clover][XTalkNum][(Crystal * (SEGS + 2)) + Seg] =
                      ((XTalkFrac[Clover][XTalkNum][(Crystal * (SEGS + 2)) + Seg] * Count) + XTalkTemp) / (Count + 1);
               }
            }
            XtalkCount[Clover][(HitCrystal * 4) + HitSeg] += 1;

         }

      }
   }

   hCloverFoldTig->Fill(CloverFoldTig);
   hCrystalFoldTig->Fill(CrystalFoldTig);
   hSegFoldTig->Fill(SegFoldTig);
   hCoreAddBackTig->Fill(CoreABTig);

}


void InitPropXtalk()
{
   char Colours[] = "BGRW";
   char name[512], title[512];
   int Clover, Crystal, Seg, Fold;

   // Initialise output file                                       
   outfile = new TFile("PropXtalkOut.root", "RECREATE");
   dRaw = outfile->mkdir("Raw");
   dSum = outfile->mkdir("Sum");
   dAddBack = outfile->mkdir("Add-Back");
   dFold = outfile->mkdir("Fold");
   dOther = outfile->mkdir("Other");
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
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         Seg = 0;
         sprintf(name, "TIG%02d%c%02dA Core En", Clover + 1, Colours[Crystal], Seg);
         sprintf(title, "TIG%02d%c%02dA Core A Energy (keV)", Clover + 1, Colours[Crystal], Seg);
         hEn[Clover][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%c%02dx Seg En", Clover + 1, Colours[Crystal], Seg);
            sprintf(title, "TIG%02d%c%02dx Seg Energy (keV)", Clover + 1, Colours[Crystal], Seg);
            hEn[Clover][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         }
         Seg = SEGS + 1;
         sprintf(name, "TIG%02d%c%02dB Core En", Clover + 1, Colours[Crystal], 0);
         sprintf(title, "TIG%02d%c%02dB Core B Energy (keV)", Clover + 1, Colours[Crystal], 0);
         hEn[Clover][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }
   dSum->cd();
   //Sum
   sprintf(name, "TIGRESS Core Sum");
   sprintf(title, "TIGRESS Core Sum Energy (keV)");
   hCoreSumTig = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      sprintf(name, "TIG%02d Core Sum", Clover + 1);
      sprintf(title, "TIG%02d Core Sum Energy (keV)", Clover + 1);
      hCoreSumClover[Clover] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      sprintf(name, "TIG%02d Seg Sum", Clover + 1);
      sprintf(title, "TIG%02d Segment Sum Energy (keV)", Clover + 1);
      hSegSumClover[Clover] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         sprintf(name, "TIG%02d%c Seg Sum", Clover + 1, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Segment Sum Energy (keV)", Clover + 1, Colours[Crystal]);
         hSegSumCrystal[Clover][Crystal] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }
   //Addback
   dAddBack->cd();
   sprintf(name, "TIGRESS Core AB");
   sprintf(title, "TIGRESS All Core Add-Back Energy (keV)");
   hCoreAddBackTig = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      sprintf(name, "TIG%02d Core AB", Clover + 1);
      sprintf(title, "TIG%02d Core Add-Back Energy (keV)", Clover + 1);
      hCoreAddBackClover[Clover] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      sprintf(name, "TIG%02d Seg AB", Clover + 1);
      sprintf(title, "TIG%02d Segment Add-Back Energy (keV)", Clover + 1);
      hSegAddBackClover[Clover] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      //Seg AB Crystal???
   }
   //Derived
   dFold->cd();
   sprintf(name, "TIG Clover Fold");
   sprintf(title, "TIGRESS Clover Fold");
   hCloverFoldTig = new TH1F(name, title, FOLD_MAX, 0, FOLD_MAX);
   sprintf(name, "TIG Crys Fold");
   sprintf(title, "TIGRESS Crystal Fold");
   hCrystalFoldTig = new TH1F(name, title, FOLD_MAX, 0, FOLD_MAX);
   sprintf(name, "TIG Seg Fold");
   sprintf(title, "TIGRESS Segment Fold");
   hSegFoldTig = new TH1F(name, title, FOLD_MAX, 0, FOLD_MAX);
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      sprintf(name, "TIG%02d Crys Fold", Clover + 1);
      sprintf(title, "TIG%02d Crystal Fold", Clover + 1);
      hCrystalFoldClover[Clover] = new TH1F(name, title, FOLD_MAX, 0, FOLD_MAX);
      sprintf(name, "TIG%02d Seg Fold", Clover + 1);
      sprintf(title, "TIG%02d Segment Fold", Clover + 1);
      hSegFoldClover[Clover] = new TH1F(name, title, FOLD_MAX, 0, FOLD_MAX);
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         sprintf(name, "TIG%02d%c Seg Fold", Clover + 1, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Segment Fold", Clover + 1, Colours[Crystal]);
         hSegFoldCrystal[Clover][Crystal] = new TH1F(name, title, FOLD_MAX, 0, FOLD_MAX);
      }
   }
   dOther->cd();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         sprintf(name, "TIG%02d%c Core En F1", Clover + 1, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Core Energy SegFold1 Any Seg (keV)", Clover + 1, Colours[Crystal]);
         hFold1CoreEn[Clover][Crystal][0] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%c Core En F1 Seg%02d", Clover + 1, Colours[Crystal], Seg);
            sprintf(title, "TIG%02d%c Core Energy SegFold1 Segment%02d (keV)", Clover + 1, Colours[Crystal], Seg);
            hFold1CoreEn[Clover][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         }
      }
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         Fold = 0;
         sprintf(name, "TIG%02d%c Seg AB", Clover + 1, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Segment Add-Back Energy Any Fold (keV)", Clover + 1, Colours[Crystal]);
         hSegAddBackCrystalByFold[Clover][Crystal][Fold] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         for (Fold = 1; Fold < SEGS; Fold++) {
            sprintf(name, "TIG%02d%c Seg AB F%01d", Clover + 1, Colours[Crystal], Fold);
            sprintf(title, "TIG%02d%c Segment Add-Back Energy Fold %01d (keV)", Clover + 1, Colours[Crystal], Fold);
            hSegAddBackCrystalByFold[Clover][Crystal][Fold] =
                new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         }
      }

      sprintf(name, "TIG%02d Seg AB F0", Clover + 1);
      sprintf(title, "TIG%02d Segment Add-Back Energy Any Fold (keV)", Clover + 1);
      hSegAddBackCloverByFold[Clover][0] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      for (Fold = 1; Fold < SEGS; Fold++) {
         sprintf(name, "TIG%02d Seg AB F%01d", Clover + 1, Fold);
         sprintf(title, "TIG%02d Segment Add-Back Energy Fold %01d (keV)", Clover + 1, Fold);
         hSegAddBackCloverByFold[Clover][Fold] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }

   // 2D crosstalk matrices 
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      sprintf(name, "TIG%02d Xtalk", Clover + 1);
      sprintf(title, "TIG%02d Fractional Crosstalk", Clover + 1);
      hXTalk[Clover] =
          new TH2F(name, title, (CRYSTALS * (SEGS + 2)), 0, (CRYSTALS * (SEGS + 2)), (CRYSTALS * (SEGS + 2)), 0,
                   (CRYSTALS * (SEGS + 2)));
   }

   // Initialise alternate gains
   if (USE_ALT_CALIB) {
      //SetGains();
   }

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
   XTalkOut.open("XTalkValuesOut.txt");
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      XTalkOut << endl << "------------------" << endl << " Clover " << Clover +
          1 << endl << "------------------" << endl << endl;
      for (HitSegment = 0; HitSegment < ((SEGS + 2) * CRYSTALS); HitSegment++) {
         for (OtherSegment = 0; OtherSegment < ((SEGS + 2) * CRYSTALS); OtherSegment++) {
            XTalkOut << XTalkFrac[Clover][HitSegment][OtherSegment] << " ";
            hXTalk[Clover]->SetBinContent(HitSegment, OtherSegment, XTalkFrac[Clover][HitSegment][OtherSegment]);
         }
         XTalkOut << endl;
      }
   }

   // Write spectra to file
   // Raw
   outfile->cd();
   dRaw->cd();
   hHitPattern->Write();
   hEHitPattern->Write();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg < SEGS + 2; Seg++) {
            hEn[Clover][Crystal][Seg]->Write();
         }
      }
   }
   // Sums
   dSum->cd();
   hCoreSumTig->Write();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      hCoreSumClover[Clover]->Write();
      hSegSumClover[Clover]->Write();
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         hSegSumCrystal[Clover][Crystal]->Write();
      }
   }

   // Add-Back
   dAddBack->cd();
   hCoreAddBackTig->Write();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      hCoreAddBackClover[Clover]->Write();
      hSegAddBackClover[Clover]->Write();
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         //hSegAddBackCrystal[Clover][Crystal]->Write();
      }
   }

   // Derived
   dFold->cd();
   hCloverFoldTig->Write();
   hCrystalFoldTig->Write();
   hSegFoldTig->Write();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      hCrystalFoldClover[Clover]->Write();      // Num crystals hit in each clover
      hSegFoldClover[Clover]->Write();
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         hSegFoldCrystal[Clover][Crystal]->Write();
      }
   }
   dOther->cd();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg < SEGS + 1; Seg++) {
            hFold1CoreEn[Clover][Crystal][Seg]->Write();
         }
      }
   }

   for (Clover = 0; Clover < CLOVERS; Clover++) {
      for (Fold = 0; Fold < SEGS; Fold++) {
         hSegAddBackCloverByFold[Clover][Fold]->Write();
      }
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Fold = 0; Fold < SEGS; Fold++) {
            hSegAddBackCrystalByFold[Clover][Crystal][Fold]->Write();
         }
      }
   }

   for (Clover = 0; Clover < CLOVERS; Clover++) {
      hXTalk[Clover]->Write();
   }

   outfile->Close();

}
