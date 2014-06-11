// Header file for stuff related to run configuration and command line options
// ---------------------------------------------------------------------------

// --------------------------------------------------------
// Definitions:
// --------------------------------------------------------

// Printing info
#define PRINT_OUTPUT 1
#define PRINT_VERBOSE 0
#define PRINT_FREQ 50000
// Main loop control
#define MAX_EVENTS 0
#define DEBUG_TREE_LOOP 0
// ROOT Stuff
#define ROOT_VIRT_SIZE    1024u*1024u*1024u     // 1x10^7 or ~10Mb seems to run fast-ish but not freeze the system completely.
                                     // that's on my 6Gb 2.6GHz i5 (YMMV)
                                     // Correction: for a longer run file 10^9 was needed
// Stuff about the experimental setup
#define CLOVERS  16
#define CRYSTALS  4
#define SEGS      8
#define INTEGRATION 125         // Integration factor applied to charge values
// Constants
#define PI 3.14159265359

//Other stuff
#define CHAR_BUFFER_SIZE 1024

// Item numbers for source records
#define SOURCE_60CO 0 
#define SOURCE_152EU_FRONT 1
#define SOURCE_152EU_BACK 2
#define SOURCE_133BA 3

// Information for calib and eff
// --------------------------------------------

// Plotting
#define PLOT_FITS 0             // plot fits, all chans plotted if below items = 0
#define PLOT_CALIB 0            // Plot calibration
#define PLOT_RESIDUAL 0
#define PLOT_CALIB_SUMMARY 1    // Plot and histo of calibration values
#define PLOT_WAVE 0

// Other info for calibration code
#define FIT_EN 1                //
#define FIT_WAVE_EN 1
#define OUTPUT_GAIN 1           // write full-run gains to file
#define OUTPUT_REPORT 1         // write full report including all fits and gains

// Definitions of types
typedef std::map < std::vector<int>, float > ReferenceValueMap;
typedef std::pair <std::vector<int>, float > ReferenceValuePair;
typedef std::map <std::vector<int>, float> :: iterator ReferenceValueMapIt;

// Define structure for storing full configuration 
struct RunConfig {              // this struct will hold all information 

   // Settings related to how the code works
   // ---------------------------------------
   // Which parts of code to run
   bool RunCalibration;
   bool RunSpecCal;
   bool RunSpecCal2;
   bool RunEfficiency;
   bool RunSpecEff;
   bool RunPropCrosstalk;
   bool RunWaveform;
   bool RunDiffCrosstalk;
   bool RunGeTiming;
   // Input files
    std::vector < std::string > files;
   // root output files
    std::string OutPath;
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
   // ROOT stuff
   unsigned int ROOT_MaxVirtSize;
   // Configuration file
   std::string ConfigFile;

   // Global physics settings
   // ------------------------------------------
   // Properties of waveforms stored in the data
   unsigned int WaveformSamples;
   unsigned int WaveInitialSamples;
   unsigned int WaveFinalSamples;
   // integration in charge evaluation
   unsigned int Integration;
   // source Information
    std::vector < std::vector < float >>Sources;
    std::vector < int >SourceNumCore;
    std::vector < int >SourceNumFront;
    std::vector < int >SourceNumBack;
   // Energy thresholds
   float EnergyThresh;
   int ChargeThresh;
   // Array mode
   bool HighEffMode;            // 1 if TIGRESS is in high eff mode (11cm), 0 if in high peak/total (14.5cm)

   // Settings for individual functions
   // ------------------------------------------

   // SortTrees:Calib() and SortHistos:CalibrateFiles() 
   // What to do
   bool CalEnergy;              // fit charge spectra
   bool CalWave;                // fit wave-charge spectra
   bool Cal2D;                  // create 2D core-seg energy matrices for low stat seg calibration
   bool CalReport;              // write full report on peak fitting as well as list of gains
   bool CalFile;                // produce a .cal file, readable by GRSISpoon
   bool CalList[CLOVERS][CRYSTALS][SEGS + 2];   // mask to determine which channels CalibSpectra() should run on
   bool CalListProvided;
   bool CalOutputBad;           // output to gain file if calibration is bad?
   bool CalOverwriteBad;           // Overwrite bad calibration values with mean (useful to get bad seg cals roughly right)
   // Output
    std::string CalOut;         // First 8 letters of this should match any file loaded in to offline cal
    std::string CalSpecOut;
    std::string CalName;
    std::string AnaName;
   bool WriteFits;              // --cal: write histo after fits, --calof: write new file with histos and fits
   // Gain drift/time dependent stuff
   float MaxTime;        // in seconds, max time after start of run.
   unsigned int TimeBins;
   float TimeBinSize;
   bool FitTempSpectra;  // should the temporary charge spectra be fitted for gain drift check
   // Charge spectra
   int ChargeBins;
   int ChargeBins2D;
   float ChargeMax;
   float WaveChargeMax;
   // What to plot
   bool PlotFits;
   bool CalibPlots[CLOVERS][CRYSTALS][SEGS + 2];        // records if fits should be plotted for each channel
   bool PlotCalib;
   bool PlotResidual;
   bool PlotCalibSummary;
   // Fit Options
   // Peak Search stuff (passed to ROOT's TSpectrum->Search() method)
   float EnSearchThresh;
   float EnSearchSigma;
   float WaveSearchThresh;
   float WaveSearchSigma;
   // Energy/ch fitting 
   // Estimate of gain for input of fit
   float EnGainEst;
   float WaveGainEst;
   // Extra calibration point at 0 ch = 0 keV
   // Peak Fitting
   unsigned int MinFitCounts;  // Minimum counts in whole spectrum for fit to be attempted
   float FitWidth_keV;           // Width of region either side of peak to be fitted
   bool FitBackground;           // 1 = yes, 0 = no.  Should be best to use this all the time but left option there just in case.
   float BackWidth_keV; 
   
   // Initial values for custom fit functions
   // These only effect the custom functions used if FIT_BACKGROUND == 1
   float EnergySigmaZero;
   float EnergySigma1MeV; 
   float WaveSigmaZero; 
   float WaveSigma1MeV;
   // Checks on fit quality
   float GausHeightMin;
   float GausCSPDMax;
   float GausSigmaMin;
   
   // Calibration options
   bool FitZero;                // Add extra calibration point at 0ch = 0keV
   float ZeroErr;
   bool ForceLinear;         // force linear calibration, even if numlines > 2
   bool ManualPeakSelect[CLOVERS][CRYSTALS][SEGS + 2];  // records if manual peak selection should be used
   bool ManualPeakCorrection;   // Manual peak selection if auto fails
   // Reference 
   ReferenceValueMap Crystal_FWHM_1332;  // # CRYSTAL_FHWM_1332

   // CoincEff()
   // What to plot
   bool PlotEff;
   bool OutputEff;
   // Output
    std::string EffOut;
    std::string EffTxtOut;
   // Reference values
   ReferenceValueMap Sim_Clover_AB_Eff;  // # SIM_CLOVER_AB_EFF
   ReferenceValueMap Exp_Clover_AB_Eff;  // # EXP_CLOVER_AB_EFF
   ReferenceValueMap Sim_Crystal_Eff;    // # SIM_CRYSTAL_EFF
   ReferenceValueMap Exp_Crystal_Eff;    // # EXP_CRYSTAL_EFF
   
   bool EffHaveSimRef;
    std::string EffSimRefFileName;
   bool EffHaveExpRef;
    std::string EffExpRefFileName;

   // PropXtalk()
   // Output
    std::string PropOut;
    std::string PropTxtOut;
    
   // GeTiming()
   // Output
   std::string GeTimingOut;
   float GeTimingGateCentre;
   float GeTimingGateWidth;
};

// --------------------------------------------------------
// Variables:
// --------------------------------------------------------

extern RunConfig Config;


// --------------------------------------------------------
// Functions:
// --------------------------------------------------------
int LoadDefaultSettings();
int ReadCommandLineSettings(int argc, char **argv);
int ReadOptionsFile(std::string filename);
int ReadCloverRef(std::ifstream *File, ReferenceValueMap *Map);
int ReadCrystalRef(std::ifstream *File, ReferenceValueMap *Map);




