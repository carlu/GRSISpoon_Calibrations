GRSISpoon_Calibrations
======================

Code to sort the ROOT tree of TTigFragments, as produced by GRSISpoon for the TIGRESS/GRIFFIN collaboration.

The idea is to provide a simple means of performing all initial calibrations for a TIGRESS or GRIFFIN experiment.  So far functions have been implemented to provide Energy, Efficiency and Crosstalk calibrations for the HPGe gamma-ray detectors, but additional functionality for ancillary detectors such as SHARC is also planned.  

To compile
----------

Same system requirements as https://github.com/GRIFFINCollaboration/GRSISpoon.
The SOURCEME script from that repo should also work to configure for compilation/running of this code.

Fragment Tree Sort:
g++ SortTrees.C CoincEff.C Calib.C PropXtalk.C CalibTools.C Options.C Utils.C -I$GRSISYS/include --std=c++0x -o SortTrees $GRSISYS/libraries/TigFormat/libFormat.so $GRSISYS/libraries/libCalManager.so $GRSISYS/libraries/libRootIOManager.so -O0 `root-config --cflags --libs` -lTreePlayer -lSpectrum -lgsl -lgslcblas -g

Histogram Sort:
To  compile: g++ SortHistos.C CalibTools.C Options.C Utils.C -I$GRSISYS/include --std=c++0x -o SortHistos $GRSISYS/libraries/TigFormat/libFormat.so $GRSISYS/libraries/libCalManager.so $GRSISYS/libraries/libRootIOManager.so -O0 `root-config --cflags --libs`  -lSpectrum -lgsl -lgslcblas -g


To run
------

./SortTrees -f InFile1 [InFile2...] [-e (energy Calibration File)] [-w (Wave calibration file)] [-s (Source)]  --cal/--prop/--eff
./SortHitos -f his??????.root/CalibOut??????.root  [-s (Source)]  --calspec/--caleff

note: Multiple files can be specified in for both types of sort.  For SortTrees all files will be summed into a signle set of spectra.  For SortHistos --calspec there should be a list of sources the same length as the list of files.  Separate fits will be done on each file assuming the sources are in the same order as the files.  A final calibration will include fits from all files.

Files and their jobs
--------------------

SortTrees.C : Builds TChain of multiple input files.  Builds index of assembled events if one is not already present in the file. Calls initialisation functions for other parts of the code.  Loops all built events passing each to any other parts of the code which active.  Calls finalisation functions.  Also contains some helper functions used elsewhere.

Calib.C : Builds charge spectra and hit patterns.  Fits core spectra at regular intervals through the run and keeps a record of results to check for gain drift.  This is a slow way to build histograms from a ROOT TTree but need to loop the events to do crosstalk stuff so may as well build energy histograms while we do so.

CoincEff.C : Performs a TIGRESS efficiency calibration using the source independent 60Co coincidence method.

PropXtalk.C : Builds calibrated energy and fold spectra.  Performs analysis of proportional crosstalk between channels in TIGRESS.

CalibTools.C : Helper functions for Calib.C.  Main caibration control functions for SortHistos.

SortHistos.C : Offline version of calib.C which carries out peak search and fits but uses the histograms output by Calib.C or the TIGRESS online analyser rather than building the spectra from scratch.  Uses the same functions from CalibTools.C so changes there should checked to work here too. 


Other Information.
------------------

CLOVERS, CRYSTALS, SEGS should refer to the number of clovers, crystals and segments in TIGRESS (or GRIFFIN) i.e. 16, 4, 8.
Clover, Crystal, Seg should refer to the number of an individual clover, crystal or segment.
Clover should always count from 1 to 16 not 0 to 15 in order to match the offical TRIUMF GRSI naming convention.
Seg may count from 0 to 9, if primary and secondary cores are included in the loop at 0 and 9.  
Crystal should also count form 1 to 4.
