GRSISpoon_Sort
==============

Code to sort the ROOT tree of TTigFragments, as produced by GRSISpoon for the TIGRESS/GRIFFIN collaboration.

Files and their jobs:

Main.C : Builds TChain of multiple input files.  Builds index of assembled events if one is not already present in the file. Calls initialisation functions for other parts of the code.  Loops all built events passing each to any other parts of the code which active.  Calls finalisation functions.  Also contains some helper functions used elesewhere.

Calib.C : Builds charge spectra and hit patterns.  Finds and identifies peaks in spectra, fits the peaks and performs a calibration.  Also fits core spectra at regular intervals through the run and keeps a record of results to check for gain drift.

CoincEff.C : Performs a TIGRESS efficiency calibration using the source independent 60Co coincidence method.

PropXtalk.C : Builds calibrated energy and fold spectra.  Performs analysis of proportional crosstalk between channels in TIGRESS.

CalibTools.C : Helper functions for Calib.C.

CalibOffline.C : Offline version of calib.C which carries out peak search and fits but uses the histograms output by Calib.C rather than building the spectra from scratch.  Uses the same functions from CalibTools.C so changes there should checked to work here too. 
