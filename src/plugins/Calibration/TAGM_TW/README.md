# TAGM timewalk calibration
Follow these steps after running `hd_root` with the **TAGM_TW** plugin.

1. `python tw.py -b <filename> <run number>`
2. `root -l -b 'display.C("results.root")'`
3. `ccdb add /PHOTON_BEAM/microscope -v default -r $RunNo-$RunNo tagm_tw_parms.out`

The macro `tw.py` will produce the timewalk correction parameters in a file called `tagm_tw_parms.out`. It will also create a ROOT file called `results.root` which can be used to check the quality of the fits. The ROOT macro 'display.C' can be used to overlay the fit and TH2I. The file `tagm_tw_parms_default.out` can be used to provide default timewalk parameters which negate the timewalk corrections.
