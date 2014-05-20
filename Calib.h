// Stuff for time based analysis
#define MAX_TIME 36000.0       // in seconds, max time after start of run.
#define TIME_BINS 20           // timed energy spectra to collect data for this many seconds
#define TIME_BIN_SIZE MAX_TIME/TIME_BINS
#define FIT_TEMP_SPECTRA 1      // 1= find and fit peaks when time bin changes.  0= do not.
#define FIT_FINAL_SPECTRA 1     // 1 = fit final spectra

// General options
#define NOTIFY_TIME_BIN 1       // print to stdout when spectra are reset
#define VERBOSE 1               // print lots of info on peak search and fitting
#define DEBUG 0                 // PRint Debugging messages to screen

// Plotting
#define PLOT_FITS 0             // plot fits, all chans plotted if below items = 0
#define PLOT_CLOVER 0           // select 1-16 to plot that only
#define PLOT_CRYSTAL 0          // 1-4
#define PLOT_SEG 0              //  1-10
#define PLOT_CALIB 0            // Plot calibration
#define PLOT_RESIDUAL 0
#define PLOT_CALIB_SUMMARY 1    // Plot and histo of calibration values

// Charge spectra stuff
//#define INTEGRATION 125 // Integration factor applied to charge values
      // This definition moved to main.h

#define MAX_LINES 12            // Maximum number of lines for any source

#define CHARGE_BINS 16384
#define CHARGE_MAX 1500000

// Use waveform energy
#define PLOT_WAVE 0
#define WAVE_CHARGE_MAX 16384

#define MAX_TOTAL_LINES 25 // maximum number of gamma lines to be fitted across all source

// Extra calibration point at 0 ch = 0 keV
#define INCLUDE_ZERO 0          // Add an extra calibration point at 0ch = 0 keV
#define FORCE_LINEAR 0          // Forces gain output file to use linear fit, even if more than 2 lines fitted
#define ZERO_ERR 0.01

// Peak Fitting
#define MIN_FIT_COUNTS 500      // Minimum counts in whole spectrum for fit to be attempted
#define FIT_WIDTH_KEV 15
#define FIT_BACKGROUND 1        // 1 = yes, 0 = no.  Should be best to use this all the time but left option there just in case.
#define BACK_WIDTH_KEV 10

// Energy/ch fitting
#define INITIAL_GAIN 0.16

struct FitSettings {
   int Source;                  // Source number
   int Integration;             // Integration used in charge evaluation
   float Dispersion;            // Dispersion of charge spectrum i.e. bins in charge spectrum / max value in charge Spectrum
   float SearchSigma;           // Sigma used by root in peak search, in charge units
   float SearchThresh;          // Threshold used by root in peak search
   float SigmaEstZero;          // used to estimate sigma for fit, in keV
   float SigmaEst1MeV;          // "                            "
   bool FitZero;                // Included 0ch=0keV in calibration?
   bool PlotOn;                 // Plot fits
   bool PeakSelect;             // Select primary peaks rather than find auto
   bool BackupPeakSelect;       // Fallback to manual peak select if auto fails.
   bool TempFit;                // Indicated this is fitting temp spectra so no plots, output, etc
   
   std::string HistName;
   std::string OutputName;
   
   
};

struct FitResult {  // Stores result of gaus+bg fit of single peak
   float Energy;
   float Const;
   float dConst;
   float Mean;
   float dMean;
   float Sigma;
   float dSigma;
   float ConstantBG;
   float dConstantBG;
   float LinearBG;
   float dLinearBG;
   float ChiSq;
   int NDF;
};

// Structure for fits for one spectrum
struct HistoFit {
   std::vector<FitResult> PeakFits;
   std::vector<bool> FitSuccess;
};

// Calibration for one spectrum
struct HistoCal {
   float LinGain[2];            // [O,G]
   float dLinGain[2];
   float LinGainFit[3];         // [O,G,CSPD]
   float dLinGainFit[2];        // [dO,dG]
   float QuadGainFit[4];        // [O,G,Q,CSPD]
   float dQuadGainFit[4];       // [dO,dG,dQ]
   int LinesUsed;
};

// Somewhere to store all the fits.

typedef std::map <float, FitResult> ChannelFitMap; // Map of FitResults for a particular channel using energy of line as key
typedef std::pair <float, FitResult> ChannelFitPair;  // PAirs with which to fill above
typedef std::map <float, FitResult>::iterator ChannelFitMapIt;  // Iterator for ChannelFitMap
typedef std::map < std::vector < int > , ChannelFitMap > MasterFitMap;  // Map to store all fits, key is vector of (Clover,Crystal,Seg)
typedef std::pair< std::vector < int > , ChannelFitMap >MasterFitPair; 
typedef std::map < std::vector < int > , ChannelFitMap >::iterator MasterFitMapIt;

// Sources
//extern float Sources[3][10];

extern TCanvas *cCalib1, *cCalib1a, *cCalib2, *ctemp;
