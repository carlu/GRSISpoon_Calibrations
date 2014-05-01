// Structure to hold TIGRESS channel name Mnemonic
struct Mnemonic {
   int arrayposition;
   int segment;
    std::string system;
    std::string subsystem;
    std::string arraysubposition;
    std::string collectedcharge;
    std::string outputsensor;
};



// Storing alternate Calibration
extern vector < string > EnCalibNames;
extern vector < vector < float >>EnCalibValues;
// Waveform calibration
extern vector < string > WaveCalibNames;
extern vector < vector < float >>WaveCalibValues;

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
