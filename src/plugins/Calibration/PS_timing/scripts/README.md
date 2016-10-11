# PS/PSC timing offset scripts

## Purpose
To determine TDC/ADC offsets and align the PS to the tagger.

## Example
To produce the offset histograms for run 11367:

`hd_root -PPLUGINS=PS_timing hd_rawdata_011367_*.evio`

The `PSCHit:DELTA_T_ADC_TDC_MAX` option should be added above and set to a large value, such as 500 ns,
if starting with offsets that are zeroed out or if the offsets are expected to change significantly.

To run the ROOT scripts:

`bash run.sh 11367 hd_root.root`

To publish the offsets of run 11367 to the ccdb for 11367 onward:

`bash publish.sh 11367 11367 inf`
