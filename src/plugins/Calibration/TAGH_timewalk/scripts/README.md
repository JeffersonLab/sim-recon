# TAGH TDC Timewalk Corrections
Follow these steps after running `hd_root` with the **TAGH_timewalk** plugin.

1. `root -b -q 'gaussian_fits.C("hd_root.root",true)'`
2. `root -b -q 'timewalk_fits.C("gaussian-fits-csv")'`
3. `ccdb add /PHOTON_BEAM/hodoscope/tdc_timewalk -v default -r $RunNo-$RunNo tdc_timewalk.txt`

The fitted histograms/graphs and distributions of the fit-parameters are saved in directories for checking that the fits make sense. After adding the timewalk parameters to the ccdb, rerun `hd_root` with the timewalk-plugin to apply the corrections to the data. Finally, repeat step 1, but with the second argument set to `false` in order to fit and print the overall timing distribution; do this for both the corrected and uncorrected histograms and compare them.
