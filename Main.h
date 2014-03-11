// --------------------------------------------------------
// Definitions:
// --------------------------------------------------------

// Which parts of code to run:
#define SORT_CALIB 1
#define SORT_OFFCAL 0
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

// Source information for calib and eff
// Peaks
#define SOURCE_NUM_FRONT 0      //1 // Source for front segments and core 0=60Co, 1=152Eu, 2=152Eu (no 121)
#define SOURCE_NUM_BACK 0       //2 // Source for back segments
#define SOURCE_NUM_CORE 0       //1 // Source for core
//#define SOURCE_NUM_FRONT 1 // Source for front segments and core 0=60Co, 1=152Eu, 2=152Eu (no 121)
//#define SOURCE_NUM_BACK 2 // Source for back segments
//#define SOURCE_NUM_CORE 1 // Source for core
#define NUM_LINES 2             // Number of lines to fit


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

struct RunConfig {  // this struct will hold all information 
                     // required to run the code

   // Input file paths
   std::vector < std::string > files;

   // Which parts of code to run
   bool RunCalibration;
   bool RunOffCal;
   bool RunEfficiency;
   bool RunPropCrosstalk;
   bool RunWaveform;
   bool RunDiffCrosstalk;
   
   // What to print and how often
   bool PrintBasic;
   int PrintFrequency;
   bool PrintVerbose;
   
   // Limit on number of events
   int EventLimit;

   // Properties of waveforms stored in the data
   int WaveformSamples;
   int WaveInitialSamples;
   int WaveFinalSamples;
   
   // Information on energy and wave calibration files
   std::string EnergyCalibrationFile;
   bool HaveAltEnergyCalibration;
   std::string WaveCalibrationFile;
   bool HaveWaveCalibration;
   
   // source Information
   std::vector< std::vector< float>> Sources;
   int SourceNumCore;
   int SourceNumFront;
   int SourceNumBack;
   
   // Plotting
   bool PlotFits;
   bool PlotCalib;
   bool PlotCalibSummary;
   
   
};

// Storing run settings
extern RunConfig Config;

// Storing alternate Calibration
extern vector<string> EnCalibNames;
extern vector<vector<float>> EnCalibValues;
// Waveform calibration
extern vector < string > WaveCalibNames;
extern vector < vector < float >> WaveCalibValues;

// tCanvas's for Calib() and CalibOffline()
extern TCanvas *cCalib1, *cCalib1a, *cCalib2, *cWave1, *ctemp;

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
