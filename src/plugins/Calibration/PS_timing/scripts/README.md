# PS/PSC Timing Offsets
Follow these steps after running `hd_root` with the **PS_timing** plugin.

The JANA configuration parameter `PSCHit:DELTA_T_ADC_TDC_MAX` should initially be set to a large value, such as 500 ns.

`run.sh` assumes that the ROOT scripts are in the same directory.

To run the ROOT scripts:

`bash run.sh 11367 hd_root.root`

To publish the results to the ccdb:

`bash publish.sh 11367`
