// Definitions related to TIGRESS efficiency measurments

#define GATE_LOW  1172
#define GATE_HIGH 1175
#define FIT_LOW 1300.0
#define FIT_HIGH 1365.0

#define MIN_FIT_COUNTS 10

#define OUTPUT_EFF 1

#define EN_SPECTRA_CHANS 8192
#define EN_SPECTRA_MAX 2048

struct FitResult {
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
