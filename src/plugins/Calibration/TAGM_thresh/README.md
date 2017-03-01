# How to determine thresholds
This plugin and set of scripts can be used to generate thresholds for the fADC, 
discriminator, and offline pulse integral cuts. In addition, it can be used to 
monitor the pulse height and integral distributions.

## Create offline threshold cuts
Simply run the plugin with the CCDB table /PHOTON_BEAM/microscope/integral_cuts 
set to 0, or with the option -PTAGMHit:CUT_FACTOR=0.

The script thresholds.py will fit the pulse integral distributions with an 
exponential and a Gaussian. The x-axis is the log of the pulse integral, which 
produces a better fit to the data when using an expo + gaus.

1. `hd_root -PPLUGINS=TAGM_thresh -PTAGMHit:CUT_FACTOR=0 /path/to/file.evio`
2. `python thresholds.py -b <root file>`

This produces the files results-thresh.root and fits-thresh.out. The first is 
a ROOT file which shows the results of the fit. The second shows a list of all 
channels with the mean pulse integral as well as the value 1-5 sigma below this 
mean.

## Making a CCDB table of threshold cuts
After running thresholds.py, a CCDB table can be made after determining what 
number of sigma would be good. 5 is typical.

1. `python MakeCCDB.py 5`

In this case, 5 sigma was used. The file tagm-thresh.txt was created which is 
in the proper format for the CCDB table /PHOTON_BEAM/microscope/integral_cuts.

## Making fADC thresholds
The CCDB table is useful for adding offline software thresholds. To make a set 
of thresholds for the fADC, generate a CCDB table. The files tagm-thresh.txt 
and the example roctagm1_spring_2017_v3.cnf files are needed.

1. `python setThresh.py`

This will generate a new fADC configuration file with thresholds based on the 
tagm-thresh.txt file. It converts from pulse integral to pulse height using a 
factor of 4.3 which was empirically measured. In addition, it adjusts to take 
into account for the baseline of 100 as well as the cable map. This is needed 
because it does not pull information from the translation table.

The resulting file is thresh.cnf.

## Making DSC thresholds
Similar to the fADC thresholds, this will create discriminator thresholds.

1. `python setThreshDSC.py`

The output file is threshDSC.cnf and uses the roctagm2_dsc_spring_2017_v3.cnf 
as a reference.
