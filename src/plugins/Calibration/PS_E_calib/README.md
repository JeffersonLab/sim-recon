# PS Energy Corrections
Follow these steps after running `hd_root` with the **PS Energy Calibration** plugin.

1. `root -l -b -q 'PSEcorr.C("hd_root.root")'`
2. `ccdb add /PHOTON_BEAM/pair_spectrometer/fine/energy_corrections -v default -r $RunNo-$RunNo Eparms-TAGM.out`

To create correction parameters run the macro with the flag CORRECTIONS set to false. The flags TAGM and TAGH can be set to true or false depending on which detector is used. It is recommended for the spring 2015 runs to use the TAGM. The output ROOT file `profiles.root` contains the histograms used to create the parameters. To verify the corrections run with the flag CORRECTIONS set to true. The histograms in check.root should be roughly horizontal if the parameters are good. For a better check, re-run the plugin with `hd_root` with the flag CORRECTIONS set to true. The initial histograms in the ROOT file should be horizontal.
