// Utility functions for GRSISpoon Calibration sorts

// GRSISpoon libraries
#include "TTigFragment.h"
#include <vector>

// --------------------------------------------------------
// Functions:
// --------------------------------------------------------

// Function to parse Mnemonic name:
void ParseMnemonic(std::string * name, Mnemonic * mnemonic);
// Convert crystal Name/number
int Col2Num(char Colour);
char Num2Col(int Crystal);
int GetDaqItemNum(int Clover,int Crystal,int Seg);
// Calibration
int ReadCalibrationFile(std::string filename, vector < string > *EnCalibNames,vector < vector < float >>*EnCalibValues);
float CalibrateEnergy(int Charge, std::vector < float >Coefficients);
// Waveform energy
float CalcWaveCharge(std::vector < int >wavebuffer);
float CalibrateWaveEnergy(float Charge, std::vector < float >Coefficients);
// Hit evaluation
int TestChargeHit(float Charge, int Integration, int Threshold);
// Save events
int SaveEvent(std::vector < TTigFragment > &ev, std::string Message);
