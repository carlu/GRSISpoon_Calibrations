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

#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "HFitInterface.h"

// TriScope libraries
//#include "TTigFragment.h"

// My libraries
#include "Options.h"
#include "SortTrees.h"
#include "Calib.h"
#include "CalibTools.h"

#include "Utils.h"

extern TApplication *App;



// FitGammaSpectrum2
// An attempt to improve low statistics fitting by fitting all peaks simulataneously 
// with gain coeffs as the parameters and peak centroids calculated from those.
int FitGammaSpectrum2(TH1F * Histo, HistoFit * Fit, HistoCal * Cal, FitSettings Settings) {
   
   // General variables 
   int Line; 
   char CharBuf[CHAR_BUFFER_SIZE];
   // Storring functions etc
   std::vector<TF1> vfPeaks;
   std::vector<ROOT::Math::WrappedMultiTF1> vwfPeaks;
   std::vector<ROOT::Fit::DataRange> vRanges;
   std::vector<ROOT::Fit::BinData> vData;
   std::vector<ROOT::Fit::Chi2Function> vChi2fPeaks;
   
   
   // Loop lines and create TF1, DataRange and BinData for each
   for(Line = 0; Line < Config.Sources.at(Settings.Source).size(); Line ++) {
   
      cout << Config.Sources.at(Settings.Source).at(Line) << " " << endl;
      
      snprintf(CharBuf,CHAR_BUFFER_SIZE,"fPeak-%d",Line);
      std::string fName = CharBuf; 
      
      // Energy = Offset + (Gain*(Charge/Int))
      // Applying gaussian formaula to peaks in spectrum, mean = measured charge / integration
      // therefore mean of gaussian = (Energy - Offset) / Gain
      // This get's complicated for quadratic calibration and also fewer variables generally 
      // makes the fit more reliable so only linear calibration here.  Fit results can
      // be reused later for a quadratic calibration
      // Previous gaus + bgd formula: "([0]*exp(-0.5*((x-[1])/[2])**2))+[3]"
      // New one:               [peak height]*exp(-0.5*((x-((Energy-[Offset])/[Gain]))/[sigma])) + [background] 
      snprintf(CharBuf,CHAR_BUFFER_SIZE,"([2]*exp(-0.5*((x- (  (%f - [0]) / [1] )  )/[3])**2))+[4]",Config.Sources.at(Settings.Source).at(Line));
      // Global fit: [0] = Offset, [1] = Gain   Unique fit: [2] = peak height [3] = sigma [4] = background
      std::string fFormula = CharBuf;
      cout << "Function Name    : " << fName << endl; 
      cout << "Function Formula : " << fFormula << endl;
      TF1 fPeak(fName.c_str(),fFormula.c_str(),0.0,Config.ChargeMax);
      vfPeaks.push_back(fPeak);
      
      ROOT::Fit::DataRange Range;
      Range.SetRange(0,Config.ChargeMax);
      
      ROOT::Fit::DataOptions opt;
      ROOT::Fit::BinData Data(opt,Range);
      ROOT::Fit::FillData(Data, Histo);
      vData.push_back(Data);
      
      // Add TF1s to WrappedMultiTF1
      ROOT::Math::WrappedMultiTF1 wfPeak(fPeak,1);
      vwfPeaks.push_back(wfPeak);
      
      // Create Chi2Functions from WrappedMultiTF1s and BinData
      ROOT::Fit::Chi2Function Chi2fPeaks(Data,wfPeak);
      vChi2fPeaks.push_back(Chi2fPeaks);
      
   }
   
   
   
   
   // Create global Chi2 function
   
   // Fitter
   
   // Set/Limit/Fix parameters
   
   
   
   return 0;
}



