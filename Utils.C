// Utility functions for GRSISpoon Calibration sorts

// C/C++ libraries:
using namespace std;
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string.h>
#include <fstream>
#include <math.h>
#include <map>

// ROOT libaries
#include <TRandom3.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TFolder.h>
#include <TH1F.h>

// GRSISpoon libraries
#include "TTigFragment.h"

// Other libraries
#include "SortTrees.h"
#include "Options.h"

// Rand for gain matching
static TRandom3 rand1;

// Function to parse Mnemonic name:
void ParseMnemonic(std::string * name, Mnemonic * mnemonic)
{
   std::string buf;
   mnemonic->system.assign(*name, 0, 2);
   mnemonic->subsystem.assign(*name, 2, 1);
   buf.clear();
   buf.assign(*name, 3, 2);
   mnemonic->arrayposition = atoi(buf.c_str());
   mnemonic->arraysubposition.assign(*name, 5, 1);
   mnemonic->collectedcharge.assign(*name, 6, 1);
   buf.clear();
   buf.assign(*name, 7, 2);
   mnemonic->segment = atoi(buf.c_str());
   mnemonic->outputsensor.assign(*name, 9, 1);
}

int Col2Num(char Colour)
{
   switch (Colour) {
   case 'B':
      return 0;
   case 'G':
      return 1;
   case 'R':
      return 2;
   case 'W':
      return 3;
   default:
      return -1;
   }
}

char Num2Col(int Crystal)
{
   switch (Crystal) {
   case 0:
      return 'B';
   case 1:
      return 'G';
   case 2:
      return 'R';
   case 3:
      return 'W';
   default:
      return 'X';
   }
}

int ReadCalibrationFile(std::string filename, vector < string > *CalibNames, vector < vector < float >>*CalibValues)
{

   if (Config.PrintBasic) {
      printf("Reading calibration file %s...\t", filename.c_str());
   }
   ifstream file;
   file.open(filename);
   if (!file) {
      if (Config.PrintBasic) {
         printf("could not open file.\n");
      }
      return -1;
   } else {
      if (Config.PrintBasic) {
         printf("File opened.\n");
      }
   }

   std::string line;
   char name[CHAR_BUFFER_SIZE], temp[CHAR_BUFFER_SIZE];
   float g0, g1, g2;
   int n = 0;

   while (getline(file, line)) {        // loop lines in file
      if (line.c_str()[0] != '#') {     // skip commented lines
         g0 = 0.0;              // reset values
         g1 = 0.0;
         g2 = 0.0;
         sscanf(line.c_str(), "%s %s %f %f %f", name, temp, &g0, &g1, &g2);     // grab name and 3 gains
         vector < float >GainTemp;
         GainTemp.push_back(g0);
         GainTemp.push_back(g1);
         GainTemp.push_back(g2);
         CalibValues->push_back(GainTemp);      // store values
         CalibNames->push_back(name);   // store name
         n += 1;                // count
      }
   }
   return n;
}


float CalibrateEnergy(int Charge, std::vector < float >Coefficients)
{

   float ChargeF = (float) Charge + rand1.Uniform();
   float TempInt = 125.0;
   float Energy = 0.0;
   if (Coefficients.size() == 0) {
      return ChargeF;
   }
   if (Config.Integration != 0) {
      TempInt = Config.Integration;
   }
   for (unsigned int i = 0; i < Coefficients.size(); i++) {
      Energy += Coefficients[i] * pow((ChargeF / TempInt), i);
   }
   return Energy;
}


float CalibrateWaveEnergy(float Charge, std::vector < float >Coefficients)
{

   Charge = Charge + rand1.Uniform();
   float Energy = 0.0;
   if (Coefficients.size() == 0) {
      return Charge;
   }
   for (unsigned int i = 0; i < Coefficients.size(); i++) {
      Energy += Coefficients[i] * pow(Charge, i);
   }
   return Energy;
}



float CalcWaveCharge(std::vector < int >wavebuffer)
{
   unsigned int Samp, Length;
   float Charge = 0.0;
   float Initial = 0.0;
   float Final = 0.0;
   Length = wavebuffer.size();

   if (wavebuffer.size() < Config.WaveInitialSamples + Config.WaveFinalSamples) {
      return 0.0;               // return 0 if wave too short
   }
   for (Samp = 0; Samp < Config.WaveInitialSamples; Samp++) {
      Initial += wavebuffer.at(Samp);
   }
   Initial /= Config.WaveInitialSamples;
   for (Samp = 0; Samp < Config.WaveFinalSamples; Samp++) {
      Final += wavebuffer.at(Length - Samp - 1);        // -1 because final sample in wbuffer seems to be spurious
   }
   Final /= Config.WaveFinalSamples;
   Charge = Final - Initial;
   return Charge;
}


// This function returns the "traditional" TIGRESS DAQ channel number for a given Cl,Cr,Seg
int GetDaqItemNum(int Clover,int Crystal,int Seg) {
   int ItemNum = -1;
   switch (Crystal) {       // Calculate channel number (old TIGRESS DAQ numbering)
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
      ItemNum = 10000;
      break;
   }
   return ItemNum;
}

int TestChargeHit(float Charge, int Integration, int Threshold) {
   if(Charge/Integration >= Threshold) {
      return 1;
   }
   else {
      return 0;
   }
}

// Save current event to root file for inspection later.
// This is useful for debugging unusual/rare events, put a gate on something
// then call this funtion to write the event out for inspection after the sort has finished.
int SaveEvent(std::vector < TTigFragment > &ev, std::string Message) {
   
   static TFile *EventOut = NULL;
   static TCanvas *cEvent = NULL;
   TDirectory *Dir = NULL;
   std::string EventOutName = Config.OutPath + Config.EventOut;
   char SpecName[CHAR_BUFFER_SIZE];
   char SpecTitle[CHAR_BUFFER_SIZE];
   char FragInfo[CHAR_BUFFER_SIZE];
   unsigned int Frag, Samp;
   TH1F *hWaveHist;
   
   // If not done already, create output file.
   if(EventOut == NULL) {
      EventOut = new TFile(EventOutName.c_str(), "RECREATE");
   }
   
   // Create folder for this event with message as title
   Dir = EventOut->mkdir(Message.c_str());
   Dir->cd();
   
   // Loop fragments in event
   for(Frag=0; Frag<ev.size(); Frag++) {
      
      // Create histogram for waveform
      // Name
      strncpy(SpecName, ev[Frag].ChannelName.c_str(),CHAR_BUFFER_SIZE);
      strncat(SpecName, " Wave",CHAR_BUFFER_SIZE);
      // Title
      strncpy(SpecTitle, ev[Frag].ChannelName.c_str(),CHAR_BUFFER_SIZE);
      strncat(SpecTitle, "Waveform ",CHAR_BUFFER_SIZE);
      // Add some other info to title
      sprintf(FragInfo, "(Q: %d, T2T: %d)", ev[Frag].Charge, ev[Frag].TimeToTrig);
      strncat(SpecTitle,FragInfo,CHAR_BUFFER_SIZE);
      // create
      hWaveHist = new TH1F(SpecName, SpecTitle, ev[Frag].wavebuffer.size(), 0, ev[Frag].wavebuffer.size());
      
      // Fill histo with wave samples
      for(Samp=0; Samp<ev[Frag].wavebuffer.size(); Samp++) {
         hWaveHist->SetBinContent(Samp,ev[Frag].wavebuffer.at(Samp));
      }
      // Write
      hWaveHist->Write();
   }
   return 0;
}



