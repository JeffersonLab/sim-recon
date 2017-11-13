# TAGM_TW usage
This plugin and scripts are designed to provide basic, initial timing and 
timewalk calibrations for the tagger microscope. Once all steps are 
complete, each TAGM channel should be aligned with the RF in both the 
TDC and ADC.

To fully align the TAGM with the rest of the detectors, use HLDetectorTiming.

**CCDB version 1.06 or greater is needed for timing.py to work.**

## Running TAGM_TW
When running `hd_root` with the **TAGM_TW** plugin, make sure to include 
the option `-PTAGMHit:DELTA_T_ADC_TDC_MAX=300` in order to pick up any 
large offsets.

`hd_root -PPLUGINS=TAGM_TW -PTAGMHit:DELTA_T_ADC_TDC_MAX=300 /path/to/file/hd_rawdata_XXXXXX_00*`

Be aware that the CCDB table /PHOTON_BEAM/microscope/integral_cuts can 
interfere with the pulse height distributions. To set these values to 0 
on the fly, run with the option -PTAGMHit:CUT_FACTOR=0.

Also, be aware that the reference RF source for this plugin is the TAGH. 
Make sure that the TAGH RF signal is calibrated before proceeding.

Follow these steps after running `hd_root`. 

## Initial ADC-RF and raw TDC-ADC calibration
This is the first step to calibrate the TAGM. Because the ADC is already 
timewalk corrected, this can be immediately aligned with an RF bucket. 
At the same time, the uncorrected TDC time can be adjusted to be aligned 
with the new ADC time. This is a rough calibration of the raw TDC time 
which will avoid the timewalk from absorbing large offsets.

**For this step, it is important to include -PTAGMHit:DELTA_T_ADC_TDC_MAX=300
when running TAGM_TW. Occasionally, the initial time difference is very large.**

1. `python timing.py -b <rootfile> <run number> rf <CCDB variation>`
2. `This produces the file adc_offsets-######.txt and tdc_offsets-######.txt`

## Timewalk corrections
The timewalk corrections are to be performed after the initial ADC-RF and 
raw TDC-ADC calibrations have been performed.

The timewalk plugin uses the TAGH RF as a reference. For a given hit, the 
RF time closest to the ADC time is selected and set as the RF time. The time 
difference between the TDC time and the RF time is taken. By doing it this 
way, there are not multiple timewalk distributions every beam period. This 
allows for large timewalks to be included in the fit.

1. `python tw.py -b <rootfile> <run number>`
2. `This creates the files results.root and tw-corr-#####.txt`
3. `root -l -b 'display.C("results.root")' to check the fits`

## Update CCDB
At this point, the TAGM timing calibrations should be finished.
The fadc, tdc, and timewalk constants can now be pushed to the CCDB.

** Note: the script push-to-ccdb.sh requires a run number WITHOUT a leading zero **
1. `./push-to-ccdb.sh <run number> <variation>
or
2. `ccdb add PHOTON_BEAM/microscope/fadc_time_offsets -v <variation> -r #-# adc_offsets-######.txt`
3. `ccdb add PHOTON_BEAM/microscope/tdc_time_offsets -v <variation> -r #-# tdc_offsets-######.txt`
4. `ccdb add PHOTON_BEAM/microscope/tdc_timewalk_corrections -v <variation> -r #-# tw-corr-######.txt`

## Check for bad CCDB values
Sometimes bad constants make their way into the CCDB. These can be found by 
running ccdbquery.py. The likely problems are wrong names for columns 101 and 102 
due to manual updates as well as unrealistically large offsets for the adc and tdc.

1. `python ccdbquery.py`

This provides 3 files: bad-cols.txt, bad-adc.txt, and bad-tdc.txt.

bad-cols.txt lists the run numbers where columns 101 and 102 are 120 and 121.

bad-adc.txt and bad-tdc.txt list the run and channel of the bad offset. 
** Note: this is channel, not column number. ** The adc/tdc tables have 122 entries.
Channel corresponds to their entry number. These are typically invidual channels as
they require more statistics to calibrate.

Using these files, the CCDB constants can be updated either by hand or by a custom script.
For the adc/tdc offsets, it is recommended to set any large offset to 7 as that is a typical
timing offset.

## Calibration validation
After all steps are complete, a calibration validation can be performed. Run 
the plugin again with the new constants.

**The constants from the previous step must be pushed to the CCDB and a new root
file created using TAGM_TW before performing this step.**

1. `python timing.py -b <rootfile> <run number> validate <CCDB variation>`
2. `This produces a file problem-channels.txt containing any errors or fits 
that are exceeding a default value. These channels should be investigated.`

## Timing resolutions
Once all of the calibrations are done, the timing resolutions can be measured. 
The timing.py script will look at the corrected TDC-RF plots and fit the 
distributions with a Gaussian. The results are provided as 1 sigma as well as 
FWHM and can be seen in the ROOT file resolutions.root.

1. `python timing.py -b <rootfile> <run number> res <CCDB variation>`

