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

// TriScope libraries
#include "TTigFragment.h"
//#include "TFSPC_Info.h"
//#include "TSharc.h"
//#include "TTigress.h"
//#include "TRf.h"
//#include "TTriFoil.h"

// My libraries
#include "CoincEff.h"
#include "Main.h"

//TStopwatch watch;

extern TApplication *App;

// Parameters and constants

// File pointers:
static TFile *outfile = 0;

static TDirectory *dEnergy, *dAddBack, *dOther = {0};

// Spectra pointers here, thinking I will start all pointer names with h
static TH1F *hTestSpectrum = 0;
static TH1F *hCrystalEn[CLOVERS][CRYSTALS] = { 0 };
static TH1F *hCloverEn[CLOVERS] = { 0 };
static TH1F *hCloverABEn[CLOVERS] = { 0 };
static TH1F *hCloverABEnGated[CLOVERS] = { 0 };
static TH1F *hCrystalEnGated[CLOVERS][CRYSTALS] = { 0 };

static TH1F *hArrayEn = 0;

// Other stuff

// Functions
void InitCoincEff();
void FinalCoincEff();
void FitPeak(TH1F * Histo, float Min, float Max, FitResult * FitRes);
//void ParseMnemonic(std::string *name,Mnemonic *mnemonic);

void CoincEff(std::vector < TTigFragment > &ev)
{
   //cout << "------New Event------- " << ev.size() << " fragments -------" << endl;
   Double_t val = 0;
   Int_t Crystal;
   Int_t Clover;
   Int_t GatePassed = 0;
   Int_t ABGatePassed = 0;
   Int_t GateCrystal = 0;
   Int_t GateClover = 0;
   Int_t ABGateClover = 0;
   Int_t i;
   float Energy;
   float CrystalEnergies[CLOVERS][CRYSTALS];
   float CloverAddBack[CLOVERS];

   memset(CrystalEnergies, 0.0, (CLOVERS * CRYSTALS * sizeof(float)));
   memset(CloverAddBack, 0.0, CLOVERS * sizeof(float));

   for (i = 0; i < ev.size(); i++) {
      Int_t slave = ((ev[i].ChannelAddress & 0x00F00000) >> 20);
      Int_t port = ((ev[i].ChannelAddress & 0x00000F00) >> 8);
      Int_t chan = (ev[i].ChannelAddress & 0x000000FF);
      std::string name = ev[i].ChannelName;

      //cout << "Slave, Port, Chan = " << slave << ", " << port << ", " << chan << "\t" << name << endl;

      Mnemonic mnemonic;
      if (name.size() >= 10) {
         ParseMnemonic(&name, &mnemonic);
      } else {
         cout << "This shouldn't happen if the odb is correctly configured!" << endl;
         continue;
      }

      // ----------------------------------------------
      // Now identify the channels and perform analysis
      // ----------------------------------------------

      // If TIGRESS
      if (mnemonic.system == "TI") {
         // Determine Crystal
         char Colour = mnemonic.arraysubposition.c_str()[0];
         Crystal = Col2Num(Colour);
         if (Crystal == -1) {
            cout << "Bad Colour: " << Colour << endl;
            continue;
         }
         // Determine Clover position
         Clover = mnemonic.arrayposition;
         //cout << "Clov,Crys = " << Clover << ", " << Crystal << endl;

         // If Core
         if (mnemonic.segment == 0 && mnemonic.outputsensor == "a") {
            //cout << "Energy: " << ev[i].ChargeCal << "\tCharge: " << ev[i].Charge << endl;
            Energy = ev[i].ChargeCal;
            CloverAddBack[Clover - 1] += Energy;
            if (Energy > Config.EnergyThresh) {
               hCloverEn[Clover - 1]->Fill(Energy);
               hCrystalEn[Clover - 1][Crystal]->Fill(Energy);
               hArrayEn->Fill(Energy);
               // Fill array with energies from this event
               CrystalEnergies[Clover - 1][Crystal] = Energy;
               // Test if gate passed
               if (Energy > GATE_LOW && Energy < GATE_HIGH) {
                  GatePassed += 1;
                  GateCrystal = Crystal;
                  GateClover = Clover;
                  hTestSpectrum->Fill(1.0);
                  //cout << endl << "Gate passed (Cl: " << Clover-1 << " Cr: " << Crystal << " En: " << CrystalEnergies[Clover-1][Crystal] << ")" << endl;
                  //cout << "Gate Clov " << GateClover << " Cr " << GateCrystal << endl;
               }
            }
            //cout << Energy << " " << GatePassed << endl;
         }
         // If segment
         else {

         }
      }
   }

   // If gate passed, loop crystals and increment gated spectra
   if (GatePassed == 1) {
      //cout << endl << "Passed!!!   Cl,Cr : ";
      for (Clover = 0; Clover < CLOVERS; Clover++) {
         for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
            //cout << Clover << ", " << Crystal << "\t";
            if (CrystalEnergies[Clover][Crystal] > Config.EnergyThresh) {
               //cout << "!!!Energy!!! ";
               if (Crystal != GateCrystal || Clover != (GateClover - 1)) {
                  //cout << "!!!Crystal!!!";
                  hCrystalEnGated[Clover][Crystal]->Fill(CrystalEnergies[Clover][Crystal]);
               }
            }
         }
      }
      //cout << endl;
   }
   // now check if gate is passed with add-back and increment gated AB spectra if so
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      if (CloverAddBack[Clover] > Config.EnergyThresh) {
         hCloverABEn[Clover]->Fill(CloverAddBack[Clover]);
         if (CloverAddBack[Clover] > GATE_LOW && CloverAddBack[Clover] < GATE_HIGH) {
            hTestSpectrum->Fill(3.0);
            ABGatePassed += 1;
            ABGateClover = Clover + 1;  // +1 to match the GateClover and GateCrystal values above
         }
      }
   }

   if (ABGatePassed == 1) {
      for (Clover = 0; Clover < CLOVERS; Clover++) {
         if ((CloverAddBack[Clover] > Config.EnergyThresh) && (Clover != (ABGateClover - 1))) {
            hCloverABEnGated[Clover]->Fill(CloverAddBack[Clover]);
         }
      }
   }

}

void InitCoincEff()
{

   char Colours[] = "BGRW";
   // Initialise output file                
   std::string tempstring = Config.OutPath + Config.EffOut;
   outfile = new TFile(tempstring.c_str(), "RECREATE");
   
   dEnergy = outfile->mkdir("Energy");
   dAddBack = outfile->mkdir("AddBack");
   dOther = outfile->mkdir("Other");

   char name[CHAR_BUFFER_SIZE], title[CHAR_BUFFER_SIZE];
   int Clover, Crystal;
   dOther->cd();
   hTestSpectrum = new TH1F("TS", "Test Spectrum", 4096, 0, 4095);
   
   dEnergy->cd();
   hArrayEn = new TH1F("TIG Sum En", "TIGRESS Sum Energy (keV)", EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      sprintf(name, "TIG%02d En", Clover + 1);
      sprintf(title, "TIG%02d Clover Energy (keV)", Clover + 1);
      hCloverEn[Clover] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         sprintf(name, "TIG%02d%c En", Clover + 1, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Core Energy (keV)", Clover + 1, Colours[Crystal]);
         if (Clover == 0 && Crystal == 0) {
            //cout << "Creating: " << name << title << endl;   
         }
         hCrystalEn[Clover][Crystal] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
         sprintf(name, "TIG%02d%c Gated", Clover + 1, Colours[Crystal]);
         sprintf(title, "TIG%02d%c Gated Core Energy (keV)", Clover + 1, Colours[Crystal]);
         hCrystalEnGated[Clover][Crystal] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      }
   }
   
   dAddBack->cd();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      sprintf(name, "TIG%02d AB", Clover + 1);
      sprintf(title, "TIG%02d Clover Add-Back Energy (keV)", Clover + 1);
      hCloverABEn[Clover] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
      sprintf(name, "TIG%02d Gated AB", Clover + 1);
      sprintf(title, "TIG%02d Gated Clover Add-Back Energy (keV)", Clover + 1);
      hCloverABEnGated[Clover] = new TH1F(name, title, EN_SPECTRA_CHANS, 0, EN_SPECTRA_MAX);
   }
}


void FinalCoincEff()
{

   int Clover, Crystal;
   float Const, dConst, Mean, dMean, Sigma, dSigma;
   float Counts, dCounts, dCountsFit, dCountsStat, Eff, dEff;
   FitResult FitRes;
   ofstream EffOut;
   char str[256];

   // Open a file to output efficiencies:   
   if (OUTPUT_EFF) {
      std::string tempstring = Config.OutPath + Config.EffTxtOut;
      EffOut.open(tempstring.c_str());
   }
   // now fit 1332.5keV peak in gain matched spectra
   memset(&FitRes, 0.0, sizeof(FitResult));
   // First the clover add-back Spectra
   if (OUTPUT_EFF) {
      EffOut << "Clover Addback: (" << hTestSpectrum->GetBinContent(4) << " counts in 1173 gate):" << endl;
   }
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      memset(&FitRes, 0.0, sizeof(FitResult));  // Clear this, should be overwritten every time but just in case...
      FitPeak(hCloverABEnGated[Clover], FIT_LOW, FIT_HIGH, &FitRes);
      Counts = (EN_SPECTRA_CHANS / EN_SPECTRA_MAX) * FitRes.Const * FitRes.Sigma * sqrt(2 * PI);
      dCountsFit = Counts * sqrt(pow(FitRes.dConst / FitRes.Const, 2) + pow(FitRes.dSigma / FitRes.Sigma, 2));  // Fitting error
      dCountsStat = sqrt(Counts);
      dCounts = Counts * sqrt(pow(dCountsFit / Counts, 2) + pow(dCountsStat / Counts, 2));
      Eff = (Counts / hTestSpectrum->GetBinContent(4)) * 100.0;
      dEff =
          Eff * sqrt(pow(dCounts / Counts, 2) +
                     pow(sqrt(hTestSpectrum->GetBinContent(4)) / hTestSpectrum->GetBinContent(4), 2));
      if (OUTPUT_EFF) {
         sprintf(str, "TIG%02d:\t", Clover + 1);
         EffOut << str << FitRes.Const << " +/- " << FitRes.dConst;
         EffOut << "\t" << FitRes.Mean << " +/- " << FitRes.dMean;
         EffOut << "\t" << FitRes.Sigma << " +/- " << FitRes.dSigma;
         EffOut << "\t" << Counts << " +/- " << dCounts;
         EffOut << "\t" << Eff << " +/- " << dEff << endl;
      }
   }
   // then crystal spectra
   if (OUTPUT_EFF) {
      EffOut << endl << "Crystals: (" << hTestSpectrum->GetBinContent(2) << " counts in 1173 gate):" << endl;
   }
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         memset(&FitRes, 0.0, sizeof(FitResult));
         FitPeak(hCrystalEnGated[Clover][Crystal], 1300.0, 1365.0, &FitRes);
         Counts = (EN_SPECTRA_CHANS / EN_SPECTRA_MAX) * FitRes.Const * FitRes.Sigma * sqrt(2 * PI);
         dCountsFit = Counts * sqrt(pow(FitRes.dConst / FitRes.Const, 2) + pow(FitRes.dSigma / FitRes.Sigma, 2));       // Fitting error
         dCountsStat = sqrt(Counts);
         dCounts = Counts * sqrt(pow(dCountsFit / Counts, 2) + pow(dCountsStat / Counts, 2));
         Eff = (Counts / hTestSpectrum->GetBinContent(4)) * 100.0;
         dEff =
             Eff * sqrt(pow(dCounts / Counts, 2) +
                        pow(sqrt(hTestSpectrum->GetBinContent(4)) / hTestSpectrum->GetBinContent(4), 2));
         if (OUTPUT_EFF) {
            sprintf(str, "TIG%02d%c:\t", Clover + 1, Num2Col(Crystal));
            EffOut << str << FitRes.Const << " +/- " << FitRes.dConst;
            EffOut << "\t" << FitRes.Mean << " +/- " << FitRes.dMean;
            EffOut << "\t" << FitRes.Sigma << " +/- " << FitRes.dSigma;
            EffOut << "\t" << Counts << " +/- " << dCounts;
            EffOut << "\t" << Eff << " +/- " << dEff << endl;
         }
      }
   }

   // Write histograms
   outfile->cd();
   dOther->cd();
   hTestSpectrum->Write();
   
   dEnergy->cd();
   hArrayEn->Write();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      hCloverEn[Clover]->Write();
      for (Crystal = 0; Crystal < CRYSTALS; Crystal++) {
         hCrystalEn[Clover][Crystal]->Write();
         hCrystalEnGated[Clover][Crystal]->Write();
      }
   }
   dAddBack->cd();
   for (Clover = 0; Clover < CLOVERS; Clover++) {
      hCloverABEn[Clover]->Write();
      hCloverABEnGated[Clover]->Write();
   }
   
   outfile->Close();
   if (OUTPUT_EFF) {
      EffOut.close();
   }
}


void FitPeak(TH1F * Histo, float Min, float Max, FitResult * FitRes)
{

   int Integral = Histo->Integral();
   //TSpectrum *Spec = new TSpectrum();
   std::string FitOptions = ("RQEM");
   int i;
   float Centre, ConstEst, BG;

   //cout << "Integral is: " << Integral << endl;
   if (Integral > MIN_FIT_COUNTS) {

      TF1 *FitRange = new TF1("GausFlatBack", "([0]*exp(-0.5*((x-[1])/[2])**2))+[3]", Min, Max);
      FitRange->SetParName(0, "Const");
      FitRange->SetParName(1, "Mean");
      FitRange->SetParName(2, "Sigma");
      FitRange->SetParName(3, "Constant Background");
      ConstEst = 0.0;
      TAxis *xaxis = Histo->GetXaxis();
      for (i = xaxis->FindBin(Min); i <= xaxis->FindBin(Max); i++) {
         if (Histo->GetBinContent(i) > ConstEst) {
            ConstEst = Histo->GetBinContent(i);
         }
      }
      FitRange->SetParameter(0, ConstEst);
      Centre = (Min + Max) / 2.0;
      FitRange->SetParameter(1, Centre);
      FitRange->SetParameter(2, ENERGY_SIGMA_ZERO + ENERGY_SIGMA_1MEV);
      BG = (Histo->GetBinContent(xaxis->FindBin(Min)) + Histo->GetBinContent(xaxis->FindBin(Min))) / 2.0;
      FitRange->SetParameter(3, BG);

      Histo->Fit(FitRange, "R,Q");

      FitRes->Const = FitRange->GetParameter(0);
      FitRes->dConst = FitRange->GetParError(0);
      FitRes->Mean = FitRange->GetParameter(1);
      FitRes->dMean = FitRange->GetParError(1);
      FitRes->Sigma = FitRange->GetParameter(2);
      FitRes->dSigma = FitRange->GetParError(2);
      FitRes->ConstantBG = FitRange->GetParameter(3);
      FitRes->dConstantBG = FitRange->GetParError(3);

      //cout << FitRes->Const << " " << FitRes->dConst << " " << FitRes->Mean << " " << FitRes->dMean << " " << FitRes->
      //  Sigma << " " << FitRes->dSigma << endl;
   }

}
