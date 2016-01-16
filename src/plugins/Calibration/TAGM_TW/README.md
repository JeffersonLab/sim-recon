# TAGM timewalk calibration
Follow these steps after running `hd_root` with the **TAGM_TW** plugin.

1. `root -l -b -q 'tw_corr.C("hd_root.root")'`
2. `ccdb add /PHOTON_BEAM/microscope -v default -r $RunNo-$RunNo tagm_tw_parms.out`

The macro `tw_corr.C` will produce the timewalk correction parameters in a file called `tagm_tw_parms.out`. It will also create a ROOT file called `results.root` which can be used to check the quality of the fits. A file called `sigmas.out` will also be generated containing the sigmas of Gaussian fits for each TAGM channel before and after the corrections. The file `tagm_tw_parms_default.out` can be used to provide default timewalk parameters which negate the timewalk corrections.
