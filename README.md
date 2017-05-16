# smalldata_tools

smalldata_tools contains code for both smalldata production at LCLS as well as code for analysis of this code. 

# smalldata production

The code for production has two main classes: SmallDataDefaultDetector and DetObject. Both are wrappers for psana classes (smalldata and the psana Detector interface respectively).  Their purpose will be described below.

SmallDataDefaultDetector:
this is a thin wrapper of the psana Detector interface. It creates two functions: one to check if a detector is in a given run (def inRun) and a common interface to the data (def getData). It is used in the smallData production code to reduce getting data of standard detectors for a given hutch to a single line that creates a list of "default" detectors. Upon start of a production job, it is checked if a given detector is present in the data (by looking at the config data) and non-present detectors are removed. For all other detectors, for each event the data is appended to a dictionnary that gets saved using the smalldata interface.
Detectors with a SmallDataDefault definition are:
lightStatus: this adds two booleans to the data which reflect the x-ray/laser on/off status. The code uses the EVR, but removed the necessity for the user to know the specific event codes used in each hutch.
ipimb box
timetool (standard fit results from DAQ, option to recalibate)
timetoolRaw (save raw time traces, option for re-fitting using ixppy derived code) 
control detector: this saved the control variables to the data and it adds the "currentStep" PV as it exists in XPP/XCS. The reason for this PV is to make a plot of somethign against a scan variable easy, even if the scan cycles though a list of positions several times (Henrik request).
acromag analog input
encoder reader
wave8
epicsPVs (using a hutch specific, but standard list)
xtcav (calls the psana xtcav code. Assume that calibrations are present. The latter are part of the "makepeds" script used in XPP/XCS)
damagedet (provides boolean for the data of a given detector is good in a given event, irrespective of the chosen fill value in smallData. This is used in the analsis code)

DetObject:
this class deals with data from detectors that are to be treated differently from experiment to experiment. At this point it includes areaDetectors and waveform Detectors, but most code deal with the areaDetector as reduction of the datasize is more crucial here.
DetObject is a different wrapper for the psanaDetector class. Upon initialization, the calibration data is saved in the object. At the beginning of the smallData production job, this data is saved to the hdf5 file so that e.g. the noise values of the used pedestal run can be used to create custom pixel masks for saved small ROIs. In the initialization step, the to be used calibration is defined. The reason for this step is that the CsPad calibration technique depends on the physics and users like to have control over that. You also add reduction algorithms to a list that will be worked through on each event. This list includes
ROI (saving sum, max, center-of-mass, optionally all pixels, addition of simple projections is also possible)
azimuthal averaging w/ geometric & polaization corrections (possibly 2-d)
droplet finding (used a 2-threshold algorithms)
photon finding 
auto-correlation of image to find speckle size (code from TJ. Slow.)
center fitting (fits the center of pixels able threshold using the user mask)

SmallDataProduction.py
main production code. User is only supposed to deal with a small subfraction of the code: at the beginning of an experiment, the DetObjects are being set up to do what will be necessary in a run. Previous users can do this themselves, new users do this with the "analysis POC" while looking at the documentation in the wiki. After the run has started, the user simply edits a function that return parameters for the data reduction as a function of the run-number (using simple else/if statements). This is done so that the run-dependent parameters are kept in the code to make reprocessing simple (e.g. you improve the calibation and would like to rerun all runs using several different ROIs because your experiment geometry changed during the course of an experiment). 
SmallDataAna_psana has functions that make getting the right parameters easy.

smallDataRun:
the submission script. Typically run in a procServ process during the experiment, but can also be used for reprocessing.
