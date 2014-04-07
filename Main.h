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
#define ROOT_VIRT_SIZE    10000000   // 1x10^7 or ~10Mb seems to run fast-ish but not freeze the system completely.
                                     // that's on my 6Gb 2.6GHz i5 (YMMV)
// Stuff about the experimental setup
#define CLOVERS  16
#define CRYSTALS  4
#define SEGS      8
#define INTEGRATION 125         // Integration factor applied to charge values
// Options for sorting
#define USE_ALT_CALIB 0
// Constants
#define PI 3.14159265359

//Other stuff
#define CHAR_BUFFER_SIZE 1024

// Information for calib and eff
// --------------------------------------------
// Peaks
#define SOURCE_NUM_FRONT 0      //1 // Source for front segments and core 0=60Co, 1=152Eu, 2=152Eu (no 121)
#define SOURCE_NUM_BACK 0       //2 // Source for back segments
#define SOURCE_NUM_CORE 0       //1 // Source for core
#define NUM_LINES 2             // Number of lines to fit

// Peak Search stuff
#define EN_SEARCH_THRESH 0.055  //0.0028  // minimum peak height for search as frac of max peak height
                                 // It seems as if this value needs to be much lower than the actual minimum peak height
#define EN_SEARCH_SIGMA 10      //50 // Expected sigma of peaks
#define WAVE_SEARCH_THRESH 0.01 //0.0028  // minimum peak height for search as frac of max peak height
                                 // It seems as if this value needs to be much lower than the actual minimum peak height
#define WAVE_SEARCH_SIGMA 20    //50 // Expected sigma of peaks
#define MIN_GAIN 0.125          // = 0.075 / 125    // These values cover from 0.5x to 2x the typical TIGRESS gain
#define MAX_GAIN 0.172          // 0.3 / 125
#define MIN_GAIN_WAVE 0.6       // = 0.075 / 125    // These values cover from 0.5x to 2x the typical TIGRESS gain
#define MAX_GAIN_WAVE 0.7       // 0.3 / 125

// Initial values for custom fit functions
// These only effect the custom functions used if FIT_BACKGROUND == 1
#define GAUS_CONST_INITIAL 100  // initial guess for peak height
#define ENERGY_SIGMA_ZERO 0.45  // sigma of peaks in keV at zero energy.
#define ENERGY_SIGMA_1MEV 0.45  // increase in sigma from zero to 1MeV, roughly same as sigma at 0, also in keV.
#define WAVE_SIGMA_ZERO 1.5     // estimated sigma, in keV, at zero, in wave emergy spectrum
#define WAVE_SIGMA_1MEV 0.4     // increase in sigma, in keV, from 0 to 1MeV
// Checks on fit quality
#define GAUS_HEIGHT_MIN 10      // minimum peak height for fit to be used in calibration
#define GAUS_CSPD_MAX 50        // maximum CSPD for peak to be used in calibration
#define GAUS_SIGMA_MIN 250      // fits with sigma below this will be ignored
                               // will use abs(sigma) for test as some peaks seem to converge on sensible but -ve sigma



// Other info for calibration code
#define FIT_EN 1                //
#define FIT_WAVE_EN 1
#define OUTPUT_GAIN 1           // write full-run gains to file
#define OUTPUT_REPORT 1         // write full report including all fits and gains


// --------------------------------------------------------
// Data:
// --------------------------------------------------------

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

struct RunConfig {              // this struct will hold all information 
   // required to run the code

   // Settings related to how the code works
   // ---------------------------------------
   // Which parts of code to run
   bool RunCalibration;
   bool RunOffCal;
   bool RunEfficiency;
   bool RunPropCrosstalk;
   bool RunWaveform;
   bool RunDiffCrosstalk;
   // Input file paths
    std::vector < std::string > files;
   // root output files
    std::string OutPath;
    std::string CalOut;         // First 8 letters of this should match any file loaded in to offline cal
    std::string CalOfOut;
    std::string EffOut;
    std::string PropOut;
   // text output files
    std::string EffTxtOut;
    std::string PropTxtOut;
   // Information on energy and wave calibration files
    std::string EnergyCalibrationFile;
   bool HaveAltEnergyCalibration;
    std::string WaveCalibrationFile;
   bool HaveWaveCalibration;
   // What to print and how often
   bool PrintBasic;
   int PrintFrequency;
   bool PrintVerbose;
   // Limit on number of events
   int EventLimit;

   // Global physics settings
   // ------------------------------------------
   // Properties of waveforms stored in the data
   int WaveformSamples;
   int WaveInitialSamples;
   int WaveFinalSamples;
   // source Information
    std::vector < std::vector < float >>Sources;
   int SourceNumCore;
   int SourceNumFront;
   int SourceNumBack;
   // Energy thresholds
   float EnergyThresh;
   int ChargeThresh;
   // Array mode
   bool HighEffMode;            // 1 if TIGRESS is in high eff mode (11cm), 0 if in high peak/total (14.5cm)


   // Physics settings for individual functions
   // ------------------------------------------
   // Calib() and CalibOffline()
   // What to do
   bool CalEnergy;              // fit charge spectra
   bool CalWave;                // fit wave-charge spectra
   bool CalReport;              // write full report on peak fitting as well as list of gains
   bool CalFile;                // produc a .cal file, readable by GRSISpoon
   // Output
   bool WriteFits;              // --cal: write histo after fits, --calof: write new file with histos and fits
   // What to plot
   bool PlotFits;
   bool PlotCalib;
   bool PlotCalibSummary;
   // CoincEff()
   // What to plot
   bool PlotEff;
   bool OutputEff;



};

// Storing run settings
extern RunConfig Config;

// Storing alternate Calibration
extern vector < string > EnCalibNames;
extern vector < vector < float >>EnCalibValues;
// Waveform calibration
extern vector < string > WaveCalibNames;
extern vector < vector < float >>WaveCalibValues;

// tCanvas's for Calib() and CalibOffline()
extern TCanvas *cCalib1, *cCalib1a, *cCalib2, *cWave1, *ctemp;

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
