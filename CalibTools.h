// Function incuded from CalibTools.C
int FitGammaSpectrum(TH1F * Histo, SpectrumFit * Fit, FitSettings Settings);
int FitSinglePeak(TH1F * Histo, int Line, float Energy, TF1 * FitRange, FitResult * FitRes, FitSettings Settings);
int CalibrationReport(SpectrumFit * Fit, ofstream & ReportOut, std::string HistName, FitSettings Settings);
