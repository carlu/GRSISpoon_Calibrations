// Utility functions for GRSISpoon Calibration sorts

// --------------------------------------------------------
// Functions:
// --------------------------------------------------------

// Function to parse Mnemonic name:
void ParseMnemonic(std::string * name, Mnemonic * mnemonic);
// Convert crystal Name/number
int Col2Num(char Colour);
char Num2Col(int Crystal);
// Calibration
float CalibrateEnergy(int Charge, std::vector < float >Coefficients);
// Waveform energy
float CalcWaveCharge(std::vector < int >wavebuffer);
float CalibrateWaveEnergy(float Charge, std::vector < float >Coefficients);
