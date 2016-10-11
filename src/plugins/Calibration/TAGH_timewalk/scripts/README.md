# TAGH TDC timewalk scripts

## Purpose
To determine TDC timewalk corrections for the tagger hodoscope.

## Example
To produce the timewalk histograms for run 11367:

`hd_root -PPLUGINS=TAGH_timewalk hd_rawdata_011367_*.evio`

To run the ROOT scripts:

`bash run.sh 11367 hd_root.root`

To publish the timewalks of run 11367 to the ccdb for run 11367:

`bash publish.sh 11367 11367 11367`
