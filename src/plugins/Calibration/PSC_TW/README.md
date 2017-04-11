# PSC_TW usage
This plugin and script are designed to provide timewalk calibrations 
for the PSC.

## tw-corr.C
This is deprecated! Use tw.py instead!

## Running PSC_TW
Follow these steps after running `hd_root` with the **PSC_TW** plugin. 

Also, be aware that the reference RF source for this plugin is the PSC. 
Make sure that the PSC RF signal is calibrated before proceeding.

## Timewalk corrections
The timewalk plugin uses the PSC RF as a reference. For a given hit, the 
RF time closest to the ADC time is selected and set as the RF time. The time 
difference between the TDC time and the RF time is taken. By doing it this 
way, there are not multiple timewalk distributions every beam period. This 
allows for large timewalks to be included in the fit.

1. `python tw.py -b <rootfile>`
2. `This creates the files results.root and tw-corr.txt`
3. `root -l -b 'display.C("results.root")' to check the fits`
4. `ccdb add /PHOTON_BEAM/pair_spectrometer/tdc_timewalk_corrections -v <variation> -r #-# tw-corr.txt`
