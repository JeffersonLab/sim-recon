#!/bin/bash

# This script pushes TAGM_TW calibrations to the CCDB

# $1 - run number (no leading 0 in the run number)
# $2 - variation

ccdb add PHOTON_BEAM/microscope/fadc_time_offsets -r $1-$1 -v $2 adc_offsets-$1.txt
ccdb add PHOTON_BEAM/microscope/tdc_time_offsets -r $1-$1 -v $2 tdc_offsets-$1.txt
ccdb add PHOTON_BEAM/microscope/tdc_timewalk_corrections -r $1-$1 -v $2 tw-corr-$1.txt
