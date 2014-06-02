// Functions for fitting a spectrum 
// --------------------------------

// Fit gamma spectrum with global gain and offset, local height and sigma, single fit function
int FitGammaSpectrumGlobal(TH1F * Histo, HistoFit * Fit, HistoCal * Cal, FitSettings Settings);

// Fit gamma spectrum with global gain and offset, local height and sigma, multi fit function with variable ranges
int FitGammaSpectrumGlobalMulti(TH1F * Histo, HistoFit * Fit, HistoCal * Cal, FitSettings Settings);


