RunNo=$1
RunNo_min=$2
RunNo_max=$3
Dir=OUTPUT/${RunNo}/offsets
ccdb add /PHOTON_BEAM/pair_spectrometer/base_time_offset -v default -r ${RunNo_min}-${RunNo_max} ${Dir}/base_time_offset.txt
ccdb add /PHOTON_BEAM/pair_spectrometer/coarse/tdc_timing_offsets -v default -r ${RunNo_min}-${RunNo_max} ${Dir}/tdc_timing_offsets_psc.txt
ccdb add /PHOTON_BEAM/pair_spectrometer/coarse/adc_timing_offsets -v default -r ${RunNo_min}-${RunNo_max} ${Dir}/adc_timing_offsets_psc.txt
ccdb add /PHOTON_BEAM/pair_spectrometer/fine/adc_timing_offsets -v default -r ${RunNo_min}-${RunNo_max} ${Dir}/adc_timing_offsets_ps.txt
