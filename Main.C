//To compile:
// g++ Main.C CoincEff.C Calib.C PropXtalk.C FitGammaSpectrum.C -I/home/cu/Code/TIGRESS/GRSISpoon/include --std=c++0x -o Sort /home/cu/Code/TIGRESS/GRSISpoon/libraries/TigFormat/libFormat.so /home/cu/Code/TIGRESS/GRSISpoon/libraries/libCalManager.so /home/cu/Code/TIGRESS/GRSISpoon/libraries/libRootIOManager.so -O2 `root-config --cflags --libs` -lTreePlayer -lSpectrum -lgsl -lgslcblas -g

// --------------------------------------------------------------------------------
// ----  Sort code for processing ROOT TTrees of TTigFragments                 ----
// ----  (as produced by GRSISPoon code from TIGRESS data)                     ----
// ----                                                                        ----
// ----  Main loop of files/trees and then events is here, sorting code        ----
// ----  included from other files.                                            ----
// ----    - C Unsworth                                                        ----
// --------------------------------------------------------------------------------


/*

Medium Term:
- Add folloowing codes for sorting:
      - CoincEff.C
      - Calib.C
      - PropXtalk.C
      - Waves.C
      - DiffXtalk.C
*/

// C/C++ libraries:
#include <iostream>
#include <unordered_set>
#include <vector>
using namespace std;
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <fstream>

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
#include <TApplication.h>
#include <TStyle.h>

// TriScope libraries
#include "TTigFragment.h"
//#include "TFSPC_Info.h"
//#include "TSharc.h"
//#include "TTigress.h"
//#include "TRf.h"
//#include "TTriFoil.h"

// My libraries
// #include "CoincEff.h"
#include "Main.h"
//#include "Gains.h"

#include "CalibrationManager.h"

// For tracking real time
TStopwatch watch;

TApplication* App;  // Pointer to root environment for plotting etc

// Counters
int FileCount  = 0;
int EventCount = 0;
int FragCount  = 0;
int BadEventCount = 0;

// Storing alternate Calibration
vector<string> CalibNames;
vector<vector<float>> CalibValues;

// Functions
void SortTree(const char* fn);
void IncSpectra();
int ReadCalibrationFile(std::string filename);

void CoincEff(std::vector<TTigFragment> &ev);
void InitCoincEff();
void FinalCoincEff();

void Calib(std::vector<TTigFragment> &ev);
void InitCalib();
void FinalCalib();

void PropXtalk(std::vector<TTigFragment> &ev);
void InitPropXtalk();
void FinalPropXtalk();

int main(int argc, char **argv) {

	gStyle->SetOptStat("iouRMen");
	
	
   //TTreePlayer *TTP = new TTreePlayer();	// This line is needed for the libraries required 
	                                       // to read the event index to be correctly loaded.
	                                       // without it the compiler thinks we don't need it and skips - CU 26 Aug 13 
	                                        
	//TFSPC_Info  *TFSPC = new TFSPC_Info(); // Similarily this causes the library for reading and storing  odb info 
	                                       // to be loaded correctly - CU 27 Aug 13        
   
   App= new TApplication("Output", 0, NULL);  // creates root environment for interacting with plots etc
   
   
   // Variables, Constants, etc
   int i,j;     
   
   // Timing
   TStopwatch StopWatch;
	StopWatch.Start();
  
   // Load any extra configuration information 
   //vector<string> CalibNames;
   //vector<vector<float>> CalibValues;
   if(USE_ALT_CALIB) {
      std::string CalFile = "Cal_run27401_w0_Quad_FixedW04.txt"; //"Cal_run27401_quad_w0.txt";
      int NumCal;
      NumCal = ReadCalibrationFile(CalFile);
      //cout << "Hello " << NumCal << endl;
      for(i=0; i<NumCal; i++) {
         cout << i ;
         cout << ": " << CalibNames.at(i); 
         cout << " " << CalibValues.at(i)[0] << " " << CalibValues[i][1] << " " << CalibValues[i][2] << endl;
      }
   }
   
   // Initialise spectra   
   if(SORT_EFF)   {
      cout << "Initialising Efficiency Spectra..." << endl;
      InitCoincEff(); 
   }
   if(SORT_CALIB) {
      cout << "Initialising Calibration Spectra..." << endl;
      InitCalib();    
   }
   if(SORT_PROP)  {
      cout << "Initialising Cross-talk Spectra..." << endl;
      InitPropXtalk();
   }   
   
   TChain *Chain = new TChain("FragmentTree");
   
   // Setup input file list                                        
   std::vector<std::string> files;
   for(i=1;i<argc;i++) {
      files.push_back(argv[i]);  
      Chain->Add(argv[i]);
      FileCount++;
   } 
   
   TTigFragment *pFrag = 0;
	std::vector<TTigFragment> evFrags;
	
	int nTrees = Chain->GetNtrees();
	int NumChainEntries = Chain->GetEntries();
	int NumChainEvents = (int)Chain->GetMaximum("TriggerId");  // This doesn't work, TrigID reset for each tree on chain.
	
	cout << "Chain Entries (frags) : " << NumChainEntries << endl;
	cout << "Chain Events          : " << NumChainEvents << endl;
	
	int TreeNum = -1;
	int LastTreeNum = -1;
	int NumTreeEntries = -1;
	int NumTreeEvents = -1;
	
	for(i=0; i<NumChainEntries; i++) {
	
	   //cout << "HERE!" << endl;
	   
	   Chain->LoadTree(i);
	   TreeNum = Chain->GetTreeNumber();
	   
	   if(TreeNum != LastTreeNum) {
	      cout << "Switching to TreeNum " << TreeNum+1 << " from " << LastTreeNum+1 << " at chain entry " << i << endl;
	      LastTreeNum = TreeNum; 
	   }
	   else {
	      continue; // This works in conjuncion with the line below "i += (NumTreeEntries - 10);"
	               // if i still falls in the same tree, then it is incremented without sorting until the next tree
	   }
	   
	    //cout << "HERE!!" << endl;
	   
	   TTree *Tree = Chain->GetTree();
	   
	   if(!Tree->GetTreeIndex())	{
		   printf("\nTree Index not found, Building index...");
		   fflush(stdout);
		   Tree->BuildIndex("TriggerId","FragmentId");
		   printf("  Done\n");
		   fflush(stdout);
      }
	   
	   TBranch *Branch = Tree->GetBranch("TTigFragment");
	   
	   //cout << "HERE!!!" << endl;
	   
	   Branch->SetAddress(&pFrag);
	   Tree->SetMaxVirtualSize(ROOT_VIRT_SIZE);
	   Branch->LoadBaskets();
	   
		NumTreeEntries = Tree->GetEntries();  
		NumTreeEvents  = (int)Tree->GetMaximum("TriggerId");		
		
		// cout << "HERE!!!! " << NumChainEntries << " " << NumTreeEntries << " " << NumTreeEvents << endl;
		
		for(int j=0;j<NumTreeEvents;j++)	{  
		   evFrags.clear();		   
		   int FragNum = 1;
		   
		   //cout << endl << "------------------------------" << endl;
		   //cout << "FragNum: " << FragNum << " j: " << j << endl;
		   
		   
		   while(Tree->GetEntryWithIndex(j,FragNum++) != -1 ) {
		      evFrags.push_back(*pFrag);
		      FragCount++;
		      if(DEBUG_TREE_LOOP) {
		         cout << "FragNum: " << FragNum << " j: " << j;
		      }		      
		   }
		   if(DEBUG_TREE_LOOP) {
		      cout << endl;
		   }
		   
		   //cout << "HERE!!!!! " << endl;
		   
		   if(FragNum < 2) {  // Init to 1, incremented on first GetEntryWithIndex, so should be at least 3 if any frags found      
		      BadEventCount++;
		      //cout << "HERE!!!!! " << BadEventCount << endl;
		      continue;
		   }
		   else {
		      EventCount++;  
		   }
		   
		   //cout << "HERE!!!!!!!!!!" << endl;
		   // do something with the evFrag vector (contains a built event)... ProcessEvent(evFrags); 
      	if(MAX_EVENTS > 0 && EventCount >= MAX_EVENTS) {break;}
         if(SORT_EFF)   {CoincEff(evFrags); } //passing vector of built events.
         if(SORT_CALIB) {Calib(evFrags);    }
         if(SORT_PROP)  {PropXtalk(evFrags);}
         if(DEBUG_TREE_LOOP) {
            cout << "ev Num =  " << j << ", ev.size() = " << evFrags.size() << endl;
         }
	    
	      // Print info to stdout
         if ( PRINT_OUTPUT && (j % PRINT_FREQ) == 0) {   
            cout << "--------------------------------------------------" << endl;  
            cout << "  Tree  " << TreeNum+1 << " / " << nTrees << endl;
            cout << "  Event " << EventCount+BadEventCount << " / " << NumChainEvents << "   (" << EventCount << " good, " << BadEventCount << " bad)" << endl;
            cout << "  Frag  " << FragCount  << " / " << NumChainEntries << endl;
            cout << "\t  in " << StopWatch.RealTime() << " seconds."  << endl;
            StopWatch.Continue();       
            cout << "--------------------------------------------------" << endl;
         }         
		   
		}
		
		if(MAX_EVENTS > 0 && EventCount >= MAX_EVENTS) {break;}
		
		Branch->DropBaskets("all");  // Clear cache before next tree	
		//i += (nEntries - 10);
		i += (NumTreeEntries - 10);  // This skips i to almost the end of the tree, any remaining entries
		                           // on this tree will be skipped by "if(TreeNum != LastTreeNum" condition above
	}
      
   // Now finalise sorts and write spectra to files....  
   if(SORT_EFF)   {FinalCoincEff();}
   if(SORT_CALIB) {FinalCalib();   }
   if(SORT_PROP)  {FinalPropXtalk();}
   
   return 0;   
}


void IncSpectra() {  // For testing only
   //hTestSpectrum->SetBinContent(1,1000);
}

// Function to parse Mnemonic name:
void ParseMnemonic(std::string *name,Mnemonic *mnemonic)	{
  std::string buf;
  mnemonic->system.assign(*name,0,2);
  mnemonic->subsystem.assign(*name,2,1);
  buf.clear();	
  buf.assign(*name,3,2);
  mnemonic->arrayposition = atoi(buf.c_str());
  mnemonic->arraysubposition.assign(*name,5,1);
  mnemonic->collectedcharge.assign(*name,6,1);
  buf.clear(); 
  buf.assign(*name,7,2);
  mnemonic->segment = atoi(buf.c_str());
  mnemonic->outputsensor.assign(*name,9,1);
}

int Col2Num(char Colour) {
   switch(Colour) {
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

char Num2Col(int Crystal) {
   switch(Crystal) {
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

int ReadCalibrationFile(std::string filename) {

   printf("Reading calibration file %s...\t", filename.c_str());
   
   ifstream file;
   file.open(filename);
   if (!file) {
      printf("could not open file.\n");
      return -1;
   }
   else {
      printf("File opened.\n");
   }
   
   
   std::string line;
   char name[16];
   float g0, g1, g2;
   int n = 0;
   
   //vector<string> CalibNames;
   //vector<vector<float>> CalibValues;
   
   while(getline(file,line)) {
      //cout << line << endl;
      if(line.c_str()[0]!='#') {
         sscanf(line.c_str(),"%s %f %f %f",name,&g0, &g1, &g2);
         //cout << name << " " << g0 << " " << g1 << " " << g2 << endl; 
         vector<float> GainTemp;
         GainTemp.push_back(g0);
         GainTemp.push_back(g1);
         GainTemp.push_back(g2);
         CalibValues.push_back(GainTemp);
         CalibNames.push_back(name);

         n += 1;
      }
   }   
   
   return n;

}

