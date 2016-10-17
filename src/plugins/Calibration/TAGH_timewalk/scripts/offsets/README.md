# TAGH timing offset scripts

## Purpose
To adjust the TDC/RF and TDC/ADC alignments prior to timewalk corrections.

## Example
To produce the offset histograms for run 11367:

`hd_root -PPLUGINS=TAGH_timewalk hd_rawdata_011367_*.evio`

To run the ROOT scripts:

`bash run.sh 11367 hd_root.root`

To publish the offsets of run 11367 to the ccdb for 11367 onward:

`bash publish.sh 11367 11367 inf`
