// Stuff for time based analysis
#define MAX_TIME 360000.0    // in seconds, max time after start of run.
#define TIME_BINS 2000  // timed energy spectra to collect data for this many seconds
#define TIME_BIN_SIZE MAX_TIME/TIME_BINS
#define FIT_TEMP_SPECTRA 0 // 1= find and fit peaks when time bin changes.  0= do not.
#define FIT_FINAL_SPECTRA 1 // 1 = fit final spectra  

// General options
#define NOTIFY_TIME_BIN 1   // print to stdout when spectra are reset
#define VERBOSE 1 // print lots of info on peak search and fitting
#define OUTPUT_GAIN 1 // write full-run gains to file
#define OUTPUT_REPORT 1 // write full report including all fits and gains
#define DEBUG 0  // PRint Debugging messages to screen

// Plotting
#define PLOT_FITS  1  // plot fits, all chans plotted if below items = 0
#define PLOT_CLOVER 6  // select 1-16 to plot that only
#define PLOT_CRYSTAL 0  // 1-4
#define PLOT_SEG 0  //  1-10
#define PLOT_CALIB 1 // Plot calibration 
#define PLOT_CALIB_SUMMARY 1 // Plot and histo of calibration values

// Charge spectra stuff
#define INTEGRATION 125 // Integration factor applied to charge values
#define CHARGE_BINS 16384  
#define CHARGE_MAX 1500000 

// Peak Search stuff
#define SEARCH_THRESH 0.0028  // minimum peak height for search as frac of max peak height
                                 // It seems as if this value needs to be much lower than the actual minimum peak height
#define SEARCH_SIGMA 50  // Expected sigma of peaks 
#define MIN_GAIN 0.125 // = 0.075 / 125    // These values cover from 0.5x to 2x the typical TIGRESS gain
#define MAX_GAIN 0.172 // 0.3 / 125  

// Peak Identification
#define SOURCE_NUM_FRONT 1 // Source for front segments and core 0=60Co, 1=152Eu, 2=152Eu (no 121)
#define SOURCE_NUM_BACK 2 // Source for back segments
#define SOURCE_NUM_CORE 1 // Source for core
#define NUM_LINES 8 // Number of lines to fit 

// Peak Fitting
#define MIN_FIT_COUNTS 500 // Minimum counts in whole spectrum for fit to be attempted
#define FIT_WIDTH 14000    // Range either side of peak to be fitted (value, not channels)
#define FIT_BACKGROUND 1   // 1 = yes, 0 = no.  Should be best to use this all the time but left option there just in case.
#define BACK_WIDTH 1000 // range for background estimate from either side of fit width. (value, not channels)

// Initial values for custom fit functions
// These only effect the custom functions used if FIT_BACKGROUND == 1
#define GAUS_CONST_INITIAL 100  // initial guess for peak height
#define GAUS_SIGMA_ZERO (2.65*INTEGRATION) // sigma of peaks at zero energy.  For Ge roughly ((1keV /2.35) / 0.15 kev/ch) * INTEGRATIOn
#define GAUS_SIGMA_1MEV (2.65*INTEGRATION) // increase in sigma from zero to 1MeV, roughly same as sigma at 0. 

// Checks on fit quality
#define GAUS_HEIGHT_MIN 10 // minimum peak height for fit to be used in calibration
#define GAUS_CSPD_MAX 10 // maximum CSPD for peak to be used in calibration
#define GAUS_SIGMA_MIN 250 // fits with sigma below this will be ignored
                               // will use abs(sigma) for test as some peaks seem to converge on sensible but -ve sigma

// Energy/ch fitting
#define INITIAL_GAIN 0.16

struct FitResult   {
   float Energy;
   float Const;
   float dConst;
   float Mean;
   float dMean; 
   float Sigma;
   float dSigma;
   float ChiSq;
   int NDF;
};

// Structure for fits and gain for one spectrum
struct SpectrumFit   {
   FitResult PeakFits[NUM_LINES];
   int FitSuccess[NUM_LINES];
   float LinGain[2];  // [O,G]
   float dLinGain[2];
   float LinGainFit[3]; // [O,G,CSPD]
   float dLinGainFit[2]; // [dO,dG]
   float QuadGainFit[4]; // [O,G,Q,CSPD]
   float dQuadGainFit[4]; // [dO,dG,dQ]
};


// Sources
extern float Sources[3][10];

extern TCanvas *c1, *c2, *ctemp;

