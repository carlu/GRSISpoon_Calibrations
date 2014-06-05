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
#include <TH3F.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

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
#include "CoincEff.h"
#include "Utils.h"


// Declare spectra
static TH1F *hEn[CLOVERS][CRYSTALS][SEGS + 2];  // energy for each individual channel, both cores and segs
static TH1F *hSegSegTime[CLOVERS][CRYSTALS]; 
static TH1F *hTimeToTrig[CLOVERS][CRYSTALS][SEGS + 2];
static TH2F *hGatedEnergyTimeToTrig;
static TH3F *hEnergyEnergyTimeToTrig;

// Other global stuff
static TFile *outfile = 0;
static TDirectory *dEnergy, *dSegSegTime, *dTimeToTrig = { 0 };


// functions
int InitGeTiming();
void FinalGeTiming();


void GeTiming(std::vector < TTigFragment > &ev) {
   
   unsigned int Frag, CalChan, Item;
   int Chan, Clover,Crystal, Seg;
   int GateClover, GateCrystal;
   int TimeDiff;
   bool GatePassed = 0;
   std::string Name;
   float En;
   char Colour;
   
   bool Hits[CLOVERS][CRYSTALS][SEGS+2] = {{{0}}};
   float Energies[CLOVERS][CRYSTALS][SEGS+2] = {{{0.0}}};
   int TimeToTrigs[CLOVERS][CRYSTALS][SEGS+2] = {{{0}}};
   int CloverSegFold[CLOVERS] = {0};
   int CrystalSegFold[CLOVERS][CRYSTALS] = {{0}};
   int CrystalFold = 0;
   int Times[2] = {0};
   float TempEn[2];
   int TempTime[2];
   bool Success = 0;
   
   //------------------------------------------------------
   // First loop fragments and store hits, times. energies
   //------------------------------------------------------
    
   for (Frag = 0; Frag < ev.size(); Frag++) {
      Chan = (ev[Frag].ChannelAddress & 0x000000FF);
      Name = ev[Frag].ChannelName;
      
      // Parse name
      Mnemonic mnemonic;
      if (Name.size() >= 10) {
         ParseMnemonic(&Name, &mnemonic);
      } else {
         cout << "Fragment Name Too Short! - This shouldn't happen if the odb is correctly configured!" << endl;
         continue;
      }
      
      // Get calibrated charge
      if (!Config.HaveAltEnergyCalibration) {
         //cout << "Using standard calibration..." << endl;
         En = ev[Frag].ChargeCal;
      } else {
         int NewCoeffFound = 0;
         std::vector < float >Coefficients;
         for (CalChan = 0; CalChan < EnCalibNames.size(); CalChan++) {
            if (strncmp(EnCalibNames[CalChan].c_str(), Name.c_str(), 9) == 0) { // bug!  this will match the first core
               // name it finds to either a OR b.  Compare 10 chars woud work but then case sensitivity isses on the x/a/b 
               // at the end.  Don't really need second core energy right now so I will come back to this later
               NewCoeffFound = 1;
               break;
            }
         }
         if (NewCoeffFound == 1) {      // If a new set of coeffs was found, then calibrate
            En = CalibrateEnergy(ev[Frag].Charge, EnCalibValues.at(CalChan));
         } else {               // else use the existing calibration
            En = ev[Frag].ChargeCal;
         }
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

         // Fill energy array and hit pattern
         if(En > Config.ChargeThresh) {
            Hits[Clover - 1][Crystal][Seg] = 1;
            hEn[Clover - 1][Crystal][Seg]->Fill(En);
            Energies[Clover - 1][Crystal][Seg] = En;
            
            TimeToTrigs[Clover - 1][Crystal][Seg] = ev[Frag].TimeToTrig;
            hTimeToTrig[Clover - 1][Crystal][Seg]->Fill(ev[Frag].TimeToTrig);
            
            // If this is a seg, count clover fold
            if(Seg > 0 && Seg < 9) {
               CloverSegFold[Clover-1] += 1;
               CrystalSegFold[Clover -1][Crystal] += 1;
            }
            if(Seg==0) {
               CrystalFold += 1;
            }
         }
      }
   }
   
   
   
   //----------------------------------------------------------------
   // Seg1->Seg2 Time difference for core photopeak segfold=2 events
   //----------------------------------------------------------------
   
   // Loop clovers
   for(Clover = 1; Clover <=CLOVERS; Clover ++) {
      Item = 0;
      // Loop crystals
      for(Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         // If energy gate passed
         if(Energies[Clover-1][Crystal][0] > (Config.GeTimingGateCentre - (Config.GeTimingGateWidth/2.0)) && 
               Energies[Clover-1][Crystal][0] < (Config.GeTimingGateCentre + (Config.GeTimingGateWidth/2.0))) {
            // If seg fold = 2
            if(CrystalSegFold[Clover -1][Crystal] == 2) {
               // Find hit segments
               for(Seg = 1; Seg <= SEGS; Seg++) {
                  if(Hits[Clover-1][Crystal][Seg] == 1) {
                     // Store TimeToTrig
                     Times[Item] = TimeToTrigs[Clover-1][Crystal][Seg];
                     Item++;
                  }
               }
               // If we found both hit segments
               if(Item==2) {
                  // Increment difference in TimeToTrigs
                  hSegSegTime[Clover-1][Crystal]->Fill(Times[0] - Times[1] + 500);
               }   
            }
         }
      }
   }
   
   // ----------------------------------------------------   
   // Gated 2D energy vs TimeToTrig
   // ----------------------------------------------------
   
   // Reset gates
   GatePassed = 0;
   GateClover = 0;
   GateCrystal = 0;
   // Loop Clovers
   for(Clover = 1; Clover <=CLOVERS; Clover ++) {
      // Loop crystals
      for(Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         // If energy gate passed
         if(Energies[Clover-1][Crystal][0] > (Config.GeTimingGateCentre - (Config.GeTimingGateWidth/2.0)) && 
            Energies[Clover-1][Crystal][0] < (Config.GeTimingGateCentre + (Config.GeTimingGateWidth/2.0))) {
            GatePassed = 1;
            GateClover = Clover;
            GateCrystal = Crystal;
         }
      }
   }
   
   // If gate passed, search for another hit crystal and record time&energy
   if(GatePassed==1) {
      // Loop clovers
      for(Clover = 1; Clover <=CLOVERS; Clover ++) {
         // Loop crystals
         for(Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            // If this is not the gate crystal
            if(Clover!=GateClover || Crystal!=GateCrystal) {
               if(Hits[Clover-1][Crystal][0] == 1) {  // If hit then this came at same time as gate
                  TimeDiff = TimeToTrigs[Clover-1][Crystal][0] - TimeToTrigs[GateClover-1][GateCrystal][0] + 500;
                  // Increment
                  hGatedEnergyTimeToTrig->Fill(Energies[Clover-1][Crystal][0], TimeDiff);
               }
            } 
         }
      }
   }
   
   // -----------------------------------------
   // 3D energy vs energy vs timetotrig diff
   // -----------------------------------------
   
   if(CrystalFold == 2) {  // Check crystalfold==2 for whole array
      // Reset values
      TempEn[0] = 0.0; TempEn[1] = 0.0;
      TempTime[0] = 0; TempTime[1] = 0;
      Item = 0; Success = 0;
      //  Loop clovers and crystals
      for(Clover = 1; Clover <=CLOVERS; Clover ++) {
         for(Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            // Check if hit, if so must be one of the two of interest
            if(Hits[Clover-1][Crystal][0] == 1) {  
               // Store energy and time
               TempEn[Item] = Energies[Clover-1][Crystal][0];
               TempTime[Item] = TimeToTrigs[Clover-1][Crystal][0];
               Item++;
            }
            // break if both ITems found
            if(Item==2) {Success = 1; break;}
         }
         // break second loop
         if(Item==2) {break;}
      }
   }
   // If success then increment 3D histogram
   if(Success == 1) {
      TimeDiff = TempTime[0]-TempTime[1] + 500;
      hEnergyEnergyTimeToTrig->Fill(TempEn[0],TempEn[1], TimeDiff);
   }
   
   
}


int InitGeTiming() {

   char name[512], title[512];
   int Clover, Crystal, Seg;
   
   // Print Gate conditions
   if(Config.PrintBasic) {
      cout << "HPGe timing sort with gate from " << (Config.GeTimingGateCentre - (Config.GeTimingGateWidth/2.0)) << " keV to ";
      cout << (Config.GeTimingGateCentre + (Config.GeTimingGateWidth/2.0)) << " keV." << endl;
   }
   
   // Initialise output file   
   std::string tempstring = Config.OutPath + Config.GeTimingOut;
   outfile = new TFile(tempstring.c_str(), "RECREATE");
   
   // Initiialise TDirectories
   dEnergy = outfile->mkdir("Energy");
   dSegSegTime = outfile->mkdir("SegSegTime");
   dTimeToTrig = outfile->mkdir("TimeToTrig");
   
   //----------------------
   // Create spectra
   //----------------------
   // Energy
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         Seg = 0;
         sprintf(name, "TIG%02d%c%02dA Core En", Clover, Num2Col(Crystal), Seg);
         sprintf(title, "TIG%02d%c%02dA Core A Energy (keV)", Clover, Num2Col(Crystal), Seg);
         hEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%c%02dx Seg En", Clover, Num2Col(Crystal), Seg);
            sprintf(title, "TIG%02d%c%02dx Seg Energy (keV)", Clover, Num2Col(Crystal), Seg);
            hEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         }
         Seg = SEGS + 1;
         sprintf(name, "TIG%02d%c%02dB Core En", Clover, Num2Col(Crystal), 0);
         sprintf(title, "TIG%02d%c%02dB Core B Energy (keV)", Clover, Num2Col(Crystal), 0);
         hEn[Clover - 1][Crystal][Seg] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }
   // Timing
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         Seg = 0;
         sprintf(name, "TIG%02d%c%02dA Time2Trig", Clover, Num2Col(Crystal), Seg);
         sprintf(title, "TIG%02d%c%02dA Time To Trigger", Clover, Num2Col(Crystal), Seg);
         hTimeToTrig[Clover - 1][Crystal][Seg] = new TH1F(name, title, 1000, 0, 1000);
         for (Seg = 1; Seg <= SEGS; Seg++) {
            sprintf(name, "TIG%02d%c%02dx Time2Trig", Clover, Num2Col(Crystal), Seg);
            sprintf(title, "TIG%02d%c%02dx Time To Trigger", Clover, Num2Col(Crystal), Seg);
            hTimeToTrig[Clover - 1][Crystal][Seg] = new TH1F(name, title, 1000, 0, 1000);
         }
         Seg = SEGS + 1;
         sprintf(name, "TIG%02d%c%02dB Time2Trig", Clover, Num2Col(Crystal), Seg);
         sprintf(title, "TIG%02d%c%02dB Time To Trigger", Clover, Num2Col(Crystal), Seg);
         hTimeToTrig[Clover - 1][Crystal][Seg] = new TH1F(name, title, 1000, 0, 1000);
      }
   }   
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         sprintf(name,"TIG%02d%c SegSeg Time",Clover, Num2Col(Crystal));
         sprintf(title, "TIG%02d%c SegSeg Time (gate = %02f keV)",Clover, Num2Col(Crystal), Config.GeTimingGateCentre);
         hSegSegTime[Clover-1][Crystal] = new TH1F(name, title, 1000, 0, 1000);
      }
   }
   
   // Matrices
   sprintf(name, "Gated En vs Time2Trig");
   sprintf(title, "Gated (%02f keV) Energy vs Time To Trigger",Config.GeTimingGateCentre);
   hGatedEnergyTimeToTrig = new TH2F(name, title, 1024, 0, EN_SPECTRA_MAX, 1000, 0, 1000);
   sprintf(name, "En vs En vs Time2Trig");
   sprintf(title, "Energy vs Energy vs Time To Trigger");
   hEnergyEnergyTimeToTrig = new TH3F(name, title, 512, 0, EN_SPECTRA_MAX, 512, 0, EN_SPECTRA_MAX, 100, 200, 800);

   return 0;
}

void FinalGeTiming() {

   int Clover, Crystal, Seg;

   // Write spectra to file
   outfile->cd();
   dEnergy->cd();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg <= SEGS + 1; Seg++) {
            hEn[Clover - 1][Crystal][Seg]->Write();
         }
      }
   }
   
   dSegSegTime->cd();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         hSegSegTime[Clover-1][Crystal] -> Write();
      }
   }
   
   dTimeToTrig->cd();
   for (Clover = 1; Clover <= CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         for (Seg = 0; Seg <= SEGS + 1; Seg++) {
            hTimeToTrig[Clover - 1][Crystal][Seg]->Write();
         }
      }
   }
   
   outfile->cd();
   hGatedEnergyTimeToTrig->Write();
   hEnergyEnergyTimeToTrig->Write();
   
   outfile->Close();
}





