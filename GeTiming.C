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

// Other global stuff
static TFile *outfile = 0;
static TDirectory *dEnergy, *dSegSegTime, *dTimeToTrig = { 0 };


// functions
int InitGeTiming();
void FinalGeTiming();


void GeTiming(std::vector < TTigFragment > &ev) {
   
   unsigned int Frag, CalChan;
   int Chan, Clover,Crystal, Seg;
   std::string Name;
   float En;
   char Colour;
   
   float Energies[CLOVERS][CRYSTALS][SEGS+2];
   int CloverSegFold[CLOVERS];
   
   
   
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
         
         hEn[Clover - 1][Crystal][Seg]->Fill(En);
         
         hTimeToTrig[Clover - 1][Crystal][Seg]->Fill(ev[Frag].TimeToTrig);
         
         // Fill energy array and hit pattern
         if(En > Config.ChargeThresh) {
            Energies[Clover - 1][Crystal][Seg] = En;
            if(Seg > 0 && Seg < 9) {
               CloverSegFold[Clover-1] += 1;
            }
         }
         
         
      }
   }
   
   // Loop clovers
   
   // Check seg fold
   
   // Check energy gate on core
   
   
   
}


int InitGeTiming() {

   char Colours[] = "BGRW";
   char name[512], title[512];
   int Clover, Crystal, Seg, Fold;
   
   // Initialise output file   
   std::string tempstring = Config.OutPath + Config.GeTimingOut;
   outfile = new TFile(tempstring.c_str(), "RECREATE");
   
   // Initiialise TDirectories
   dEnergy = outfile->mkdir("Energy");
   dSegSegTime = outfile->mkdir("SegSegTime");
   dTimeToTrig = outfile->mkdir("TimeToTrig");
   
   // Create spectra
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
         sprintf(title, "TIG%02d%c SegSeg Time",Clover, Num2Col(Crystal));
         hSegSegTime[Clover-1][Crystal] = new TH1F(name, title, 1000, 0, 1000);
      }
   }

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
   outfile->Close();
}





