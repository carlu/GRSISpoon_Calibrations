// --------------------------------------------------------
// Definitions:
// --------------------------------------------------------

// Which parts of code to run:
#define SORT_CALIB 1
#define SORT_EFF   0
#define SORT_WAVES 0
#define SORT_PROP  0
#define SORT_DIFF  0
// Printing info
#define PRINT_OUTPUT 1
#define PRINT_VERBOSE 0
#define PRINT_FREQ 50000
// Main loop control
#define MAX_EVENTS 0
#define DEBUG_TREE_LOOP 0
// ROOT Stuff
#define ROOT_VIRT_SIZE    500000000  //  500MB 
// Stuff about the experimental setup
#define CLOVERS  16
#define CRYSTALS  4
#define SEGS      8
#define INTEGRATION 125 // Integration factor applied to charge values
// Options for sorting
#define EN_THRESH 2  // energies less than this keV ignored
#define USE_ALT_CALIB 0
// Constants
#define PI 3.14159265359

// Calculation of Energy from waveform
#define WAVE_SAMPS 200
#define INITIAL_SAMPS 65
#define FINAL_SAMPS 65

//Other stuff
#define CHAR_BUFFER_SIZE 1024

// --------------------------------------------------------
// Data:
// --------------------------------------------------------

// Structure to hold TIGRESS channel name Mnemonic
struct Mnemonic	{
   int arrayposition;
   int	segment;
   std::string system;
   std::string subsystem;
   std::string arraysubposition;
   std::string collectedcharge;
   std::string outputsensor;
};

struct RunConfig {

   std::vector < std::string > files;

   bool RunCalibration;
   bool RunEfficiency;
   bool RunPropCrosstalk;
   bool RunWaveform;
   bool RunDiffCrosstalk;
   
   bool PrintBasic;
   int PrintFrequency;
   bool PrintVerbose;
   
   int EventLimit;

   int WaveformSamples;
   int WaveInitialSamples;
   int WaveFinalSamples;
   
   std::string EnergyCalibrationFile;
   bool UseAltEnergyCalibration;
   std::string WaveCalibrationFile;
   bool HaveWaveCalibration;
   
   std::vector< std::vector< float>> Sources;
};

// Storing run settings
extern RunConfig Config;

// Storing alternate Calibration
extern vector<string> EnCalibNames;
extern vector<vector<float>> EnCalibValues;
// Waveform calibration
extern vector < string > WaveCalibNames;
extern vector < vector < float >> WaveCalibValues;

//extern TCanvas *cWave1, *ctemp;

// --------------------------------------------------------
// Functions:
// --------------------------------------------------------

// Function to parse Mnemonic name:
void ParseMnemonic(std::string *name,Mnemonic *mnemonic);
// Convert crystal Name/number
int Col2Num(char Colour);
char Num2Col(int Crystal);
// Calibration
float CalibrateEnergy(int Charge, std::vector<float> Coefficients);
// Waveform energy
float CalcWaveCharge(std::vector<int>  wavebuffer);
float CalibrateWaveEnergy(float Charge, std::vector<float> Coefficients);
