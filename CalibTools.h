


//Functions
int CalibrateFiles(); 
int FitHistoFile(TFile *file, int FileType, int FileNum, MasterFitMap *FitMap, MasterFitMap *WaveFitMap);
int CalibrateChannel(ChannelFitMap Fits, FitSettings Settings);
int ConfigureEnergyFit(int Clover, int Crystal, int Seg, int FileType, int FileNum, FitSettings *Settings) ;
int ConfigureWaveEnFit(int Clover, int Crystal, int Seg,  int FileType, int FileNum, FitSettings *Settings) ;

// Function incuded from CalibTools.C
int FitGammaSpectrum(TH1F * Histo, HistoFit * Fit, FitSettings Settings);
int FitSinglePeak(TH1F * Histo, int Line, float Energy, TF1 * FitRange, FitResult * FitRes, FitSettings Settings);
int CalibrationReport(HistoFit * Fit, ofstream & ReportOut, std::string HistName, FitSettings Settings);
int WriteCalFile(HistoFit *Fit, ofstream &CalFileOut, std::string HistName, FitSettings Settings);
