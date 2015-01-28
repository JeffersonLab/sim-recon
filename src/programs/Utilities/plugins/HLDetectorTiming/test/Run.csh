#!/bin/csh

# Set up the environment you will use for the running
source /home/mstaib/commissioning/sim-recon-commissioning/Linux_RHEL6-x86_64-gcc4.4.7/setenv.csh
setenv HALLD_MY /home/mstaib/gluex/

# Set the run number and directory where to find the files
#set run=2120
set run=2209
#set directory=/raid12/gluex/fcal_bcal_m8/2trackskim/
set directory=/raid12/gluex/rawdata/Run002209/

# This will find all of the evio files in the directory matching the run number
# Don't ask me why this is how to sort it @_@ This way the run tags appear sequentially
set allfiles=`find ${directory} -maxdepth 1 -name "*${run}*.evio" |xargs ls `

# We want to attempt to do the calibration on a SQLite version of the CCDB -nc will force no overwrite of the current CCDB
wget --no-check-certificate -nc https://halldweb1.jlab.org/dist/ccdb.sqlite
setenv CCDB_CONNECTION sqlite:///`pwd`/ccdb.sqlite
setenv JANA_CALIB_URL sqlite:///`pwd`/ccdb.sqlite

# Run the plugin to see where the CCDB calibrations are before calibrations
set nThreads=32
# Measure Baseline
#set nEvents=500000
#hd_root -PPLUGINS=HLDetectorTiming,monitoring_hists -PBFIELD_MAP=Magnets/Solenoid/solenoid_1200A_poisson_20140520 $allfiles -PNTHREADS=$nThreads -PEVENTS_TO_KEEP=$nEvents -PHLDETECTORTIMING:DO_VERIFY=1
#mv hd_root.root BeforeCalibration.root

# Now we need to zero out the entries we will be calibrating so we don't have random starting points
echo "Zeroing out all relevant CCDB tables"

echo "0.0" > oneVal.txt
echo "0.0 0.0" > twoVal.txt

ccdb add CDC/base_time_offset -r ${run}-${run} oneVal.txt
ccdb add FCAL/base_time_offset -r ${run}-${run} oneVal.txt
ccdb add FDC/base_time_offset -r ${run}-${run} twoVal.txt
ccdb add START_COUNTER/base_time_offset -r ${run}-${run} twoVal.txt
ccdb add TOF/base_time_offset -r ${run}-${run} twoVal.txt
ccdb add BCAL/base_time_offset -r ${run}-${run} twoVal.txt
ccdb add PHOTON_BEAM/microscope/base_time_offset -r ${run}-${run} twoVal.txt
ccdb add PHOTON_BEAM/hodoscope/base_time_offset -r ${run}-${run} twoVal.txt

# There are some more complicated tables that need to be zero'd out as well. 
# This is slightly more complicated since each detector has it's own scheme for how the constants are stored.
# It can be done! (in hopefully less than 3000 lines)

# CDC
set CDCCounter=1
echo "0.0" >! cdc_timing_offsets.txt
while ($CDCCounter < 3522)
    echo "0.0" >> cdc_timing_offsets.txt
    @ CDCCounter++
end
ccdb add CDC/timing_offsets -r ${run}-${run} cdc_timing_offsets.txt

# FDC
# Skip this one for right now, have to understand the constants better

# SC ADC/TDC
set SCCounter=1
echo "0.0" >! sc_tdc_timing_offsets.txt
echo "0.0" >! sc_adc_timing_offsets.txt
while ($SCCounter < 30)
    echo "0.0" >> sc_tdc_timing_offsets.txt
    echo "0.0" >> sc_adc_timing_offsets.txt
    @ SCCounter++
end
ccdb add START_COUNTER/tdc_timing_offsets -r ${run}-${run} sc_tdc_timing_offsets.txt
ccdb add START_COUNTER/adc_timing_offsets -r ${run}-${run} sc_adc_timing_offsets.txt

# BCAL ADC
set BCALCounter=1
echo "0.0" >! bcal_adc_timing_offsets.txt
while ($BCALCounter < 1536)
    echo "0.0" >> bcal_adc_timing_offsets.txt
    @ BCALCounter++
end
ccdb add BCAL/ADC_timing_offsets -r ${run}-${run} bcal_adc_timing_offsets.txt

#TOF ADC/TDC
set TOFCounter=1
echo "0.0" >! tof_tdc_timing_offsets.txt
echo "0.0" >! tof_adc_timing_offsets.txt
while ($TOFCounter < 176)
    echo "0.0" >> tof_tdc_timing_offsets.txt
    echo "0.0" >> tof_adc_timing_offsets.txt
    @ TOFCounter++
end
ccdb add TOF/timing_offsets -r ${run}-${run} tof_tdc_timing_offsets.txt
ccdb add TOF/adc_timing_offsets -r ${run}-${run} tof_adc_timing_offsets.txt

# FCAL ADC
set FCALCounter=1
echo "0.0" >! fcal_adc_timing_offsets.txt
while ($FCALCounter < 2800)
    echo "0.0" >> fcal_adc_timing_offsets.txt
    @ FCALCounter++
end
ccdb add FCAL/timing_offsets -r ${run}-${run} fcal_adc_timing_offsets.txt

# TAGH ADC/TDC
# Here the data has two columns, so things are a bit more interesting
set TAGHCounter=2
echo "1 0.0" >! tagh_tdc_timing_offsets.txt
echo "1 0.0" >! tagh_adc_timing_offsets.txt
while ($TAGHCounter <= 274)
    echo "$TAGHCounter 0.0" >> tagh_tdc_timing_offsets.txt
    echo "$TAGHCounter 0.0" >> tagh_adc_timing_offsets.txt
    @ TAGHCounter++
end
ccdb add PHOTON_BEAM/hodoscope/tdc_time_offsets -r ${run}-${run} tagh_tdc_timing_offsets.txt
ccdb add PHOTON_BEAM/hodoscope/fadc_time_offsets -r ${run}-${run} tagh_adc_timing_offsets.txt

# TAGM ADC/TDC
# This one is the most fun, 3 columns
set TAGMColumnCounter=2
echo "0 1 0.0" >! tagm_tdc_timing_offsets.txt
echo "0 1 0.0" >! tagm_adc_timing_offsets.txt
while ($TAGMColumnCounter <= 102)
    echo "0 $TAGMColumnCounter 0.0" >> tagm_tdc_timing_offsets.txt
    echo "0 $TAGMColumnCounter 0.0" >> tagm_adc_timing_offsets.txt
    set TAGMRowCounter=1
    if ( $TAGMColumnCounter == 7 || $TAGMColumnCounter == 25 || $TAGMColumnCounter == 79 || $TAGMColumnCounter == 97 ) then
        while ($TAGMRowCounter <= 5)
            echo "$TAGMRowCounter $TAGMColumnCounter 0.0" >> tagm_tdc_timing_offsets.txt
            echo "$TAGMRowCounter $TAGMColumnCounter 0.0" >> tagm_adc_timing_offsets.txt
            @ TAGMRowCounter++
        end
    endif
    @ TAGMColumnCounter++
end
ccdb add PHOTON_BEAM/microscope/tdc_time_offsets -r ${run}-${run} tagm_tdc_timing_offsets.txt
ccdb add PHOTON_BEAM/microscope/fadc_time_offsets -r ${run}-${run} tagm_adc_timing_offsets.txt

rm oneVal.txt
rm twoVal.txt

# We can make our first pass at the data, just making rough timing corrections
set nEvents=100000

hd_root -PPLUGINS=HLDetectorTiming -PBFIELD_MAP=Magnets/Solenoid/solenoid_1200A_poisson_20140520 $allfiles -PNTHREADS=$nThreads -PEVENTS_TO_KEEP=$nEvents -PHLDETECTORTIMING:DO_ROUGH_TIMING=1
mv hd_root.root RoughTiming.root

# Now take the constants that were calculated and add them to the CCDB

ccdb add CDC/base_time_offset -r ${run}-${run} cdc_base_time.txt
ccdb add FCAL/base_time_offset -r ${run}-${run} fcal_base_time.txt
ccdb add FDC/base_time_offset -r ${run}-${run} fdc_base_time.txt
ccdb add START_COUNTER/base_time_offset -r ${run}-${run} sc_base_time.txt
ccdb add TOF/base_time_offset -r ${run}-${run} tof_base_time.txt
ccdb add BCAL/base_time_offset -r ${run}-${run} bcal_base_time.txt
ccdb add PHOTON_BEAM/microscope/base_time_offset -r ${run}-${run} tagm_base_time.txt
ccdb add PHOTON_BEAM/hodoscope/base_time_offset -r ${run}-${run} tagh_base_time.txt

# Next run the TDC/ADC alignment with these results
set nEvents=250000

hd_root -PPLUGINS=HLDetectorTiming -PBFIELD_MAP=Magnets/Solenoid/solenoid_1200A_poisson_20140520 $allfiles -PNTHREADS=$nThreads -PEVENTS_TO_KEEP=$nEvents -PHLDETECTORTIMING:DO_TDC_ADC_ALIGN=1
mv hd_root.root TDCADCAlign.root

ccdb add START_COUNTER/tdc_timing_offsets -r ${run}-${run} sc_tdc_timing_offsets.txt
ccdb add TOF/timing_offsets -r ${run}-${run} tof_tdc_timing_offsets.txt
ccdb add PHOTON_BEAM/hodoscope/tdc_time_offsets -r ${run}-${run} tagh_tdc_timing_offsets.txt
ccdb add PHOTON_BEAM/microscope/tdc_time_offsets -r ${run}-${run} tagm_tdc_timing_offsets.txt

exit

# Do the Track Based alignment
set nEvents=500000
hd_root -PPLUGINS=HLDetectorTiming -PBFIELD_MAP=Magnets/Solenoid/solenoid_1200A_poisson_20140520 $allfiles -PNTHREADS=$nThreads -PEVENTS_TO_KEEP=$nEvents -PHLDETECTORTIMING:DO_TRACK_BASED=1 -PTRKFIT:MASS_HYPOTHESES_POSITIVE=0.14 -PTRKFIT:MASS_HYPOTHESES_NEGATIVE=0.14 
mv hd_root.root TrackBased.root

ccdb add PHOTON_BEAM/hodoscope/tdc_time_offsets -r ${run}-${run} tagh_tdc_timing_offsets.txt
ccdb add PHOTON_BEAM/hodoscope/fadc_time_offsets -r ${run}-${run} tagh_adc_timing_offsets.txt
ccdb add PHOTON_BEAM/microscope/tdc_time_offsets -r ${run}-${run} tagm_tdc_timing_offsets.txt
ccdb add PHOTON_BEAM/microscope/fadc_time_offsets -r ${run}-${run} tagm_adc_timing_offsets.txt
ccdb add FCAL/base_time_offset -r ${run}-${run} fcal_base_time.txt
ccdb add TOF/base_time_offset -r ${run}-${run} tof_base_time.txt
ccdb add BCAL/base_time_offset -r ${run}-${run} bcal_base_time.txt

# Verify Results
set nEvents=500000
hd_root -PPLUGINS=HLDetectorTiming,monitoring_hists -PBFIELD_MAP=Magnets/Solenoid/solenoid_1200A_poisson_20140520 $allfiles -PNTHREADS=$nThreads -PEVENTS_TO_KEEP=$nEvents -PHLDETECTORTIMING:DO_VERIFY=1 
mv hd_root.root FinalResult.root

