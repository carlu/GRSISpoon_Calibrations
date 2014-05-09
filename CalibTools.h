// Function incuded from CalibTools.C
int FitGammaSpectrum(TH1F * Histo, HistoFit * Fit, FitSettings Settings);
int FitSinglePeak(TH1F * Histo, int Line, float Energy, TF1 * FitRange, FitResult * FitRes, FitSettings Settings);
int CalibrateGammaSpectrum(HistoFit * Fit, FitSettings Settings);
int CalibrationReport(HistoFit * Fit, ofstream & ReportOut, std::string HistName, FitSettings Settings);
int WriteCalFile(HistoFit *Fit, ofstream &CalFileOut, std::string HistName, FitSettings Settings);
