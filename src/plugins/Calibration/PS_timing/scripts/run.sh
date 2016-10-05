RunNo=$1
OutputDir=OUTPUT/${RunNo}
InputFile=$2
cp -p ${InputFile} hd_root.root
mkdir -p offsets
ccdb dump /PHOTON_BEAM/pair_spectrometer/base_time_offset:${RunNo} > offsets/base_time_offset_ccdb.txt
ccdb dump /PHOTON_BEAM/pair_spectrometer/coarse/tdc_timing_offsets:${RunNo} > offsets/tdc_timing_offsets_psc_ccdb.txt
ccdb dump /PHOTON_BEAM/pair_spectrometer/coarse/adc_timing_offsets:${RunNo} > offsets/adc_timing_offsets_psc_ccdb.txt
ccdb dump /PHOTON_BEAM/pair_spectrometer/fine/adc_timing_offsets:${RunNo} > offsets/adc_timing_offsets_ps_ccdb.txt
root -b -q 'fits.C("hd_root.root",true)'
root -b -q 'offsets.C("fits-csv")'
rm -f hd_root.root
mkdir -p ${OutputDir}
mv offsets ${OutputDir}; mv fits-csv ${OutputDir}; mv fits_* ${OutputDir}
