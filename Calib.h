// General options
#define DEBUG 0                 // PRint Debugging messages to screen

// Charge spectra stuff
//#define INTEGRATION 125 // Integration factor applied to charge values
      // This definition moved to main.h

#define MAX_LINES 12            // Maximum number of lines for any source
#define MAX_TOTAL_LINES 25      // maximum number of gamma lines to be fitted across all source

struct FitSettings {
   int Source;                  // Source number
   int Integration;             // Integration used in charge evaluation
   float Dispersion;            // Dispersion of charge spectrum i.e. bins in charge spectrum / max value in charge Spectrum
   float SearchSigma;           // Sigma used by root in peak search, in charge units
   float SearchThresh;          // Threshold used by root in peak search
   float SigmaEstZero;          // used to estimate sigma for fit, in keV
   float SigmaEst1MeV;          //           "          
   float GainEst;               //           "
   bool FitZero;                // Included 0ch=0keV in calibration?
   bool PlotOn;                 // Plot fits
   bool PeakSelect;             // Select primary peaks rather than find auto
   bool BackupPeakSelect;       // Fallback to manual peak select if auto fails.
   bool TempFit;                // Indicated this is fitting temp spectra so no plots, output, etc

    std::string HistName;
    std::string OutputName;


};

struct FitResult {              // Stores result of gaus+bg fit of single peak
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
   std::vector < FitResult > PeakFits;
   std::vector < bool > FitSuccess;
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

typedef std::map < float, FitResult > ChannelFitMap;    // Map of FitResults for a particular channel using energy of line as key
typedef std::pair < float, FitResult > ChannelFitPair;  // PAirs with which to fill above
typedef std::map < float, FitResult >::iterator ChannelFitMapIt;        // Iterator for ChannelFitMap
typedef std::map < std::vector < int >, ChannelFitMap > MasterFitMap;   // Map to store all fits, key is vector of (Clover,Crystal,Seg)
typedef std::pair < std::vector < int >, ChannelFitMap > MasterFitPair;
typedef std::map < std::vector < int >, ChannelFitMap >::iterator MasterFitMapIt;

// Sources
//extern float Sources[3][10];

extern TCanvas *cCalib1, *cCalib1a, *cCalib2, *ctemp;
