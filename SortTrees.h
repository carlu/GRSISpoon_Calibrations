// Structure to hold TIGRESS channel name Mnemonic
struct Mnemonic {
   int arrayposition;
   int segment;
    std::string system;
    std::string subsystem;
    std::string arraysubposition;
    std::string collectedcharge;
    std::string outputsensor;
};


// Storing alternate Calibration
extern vector < string > EnCalibNames;
extern vector < vector < float >>EnCalibValues;
// Waveform calibration
extern vector < string > WaveCalibNames;
extern vector < vector < float >>WaveCalibValues;


