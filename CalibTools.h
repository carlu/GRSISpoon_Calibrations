// Functions for managing calibration in SortHistos
// ------------------------------------------------
// Loop all files/sources
int CalibrateFiles();
// Loop histos in one file
int FitHistoFile(TFile * file, int FileType, int FileNum, MasterFitMap * FitMap, MasterFitMap * WaveFitMap);

// Functions for fitting a spectrum 
// --------------------------------
// Find/id first two peaks, rough calibration, loop all peaks
int FitGammaSpectrum(TH1F * Histo, HistoFit * Fit, HistoCal * Cal, FitSettings Settings);
// Fit a peak
int FitSinglePeak(TH1F * Histo, int Line, float Energy, TF1 * FitRange, FitResult * FitRes, FitSettings Settings);

// Functions for calibration
// --------------------------
// Calibrate a channel
int CalibrateChannel(ChannelFitMap Fits, FitSettings Settings, HistoFit * Fit, HistoCal * Cal);

// Other general helper functions
// ------------------------------
// Configure fit settings for Energy/WaveEn spectrum
int ConfigureEnergyFit(int Clover, int Crystal, int Seg, int FileType, int FileNum, FitSettings * Settings);
int ConfigureWaveEnFit(int Clover, int Crystal, int Seg, int FileType, int FileNum, FitSettings * Settings);
// Generate text report of fits and calibrations
int CalibrationReport(HistoFit * Fit, HistoCal * Cal, ofstream & ReportOut, std::string HistName, FitSettings Settings);
// Write out .cal file (not implemented yet)
int WriteCalFile(HistoFit * Fit, ofstream & CalFileOut, std::string HistName, FitSettings Settings);
