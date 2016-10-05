RunNo=$1
Dir=OUTPUT/${RunNo}/offsets
ccdb add /PHOTON_BEAM/pair_spectrometer/base_time_offset -v default -r ${RunNo}-${RunNo} ${Dir}/base_time_offset.txt
ccdb add /PHOTON_BEAM/pair_spectrometer/coarse/tdc_timing_offsets -v default -r ${RunNo}-${RunNo} ${Dir}/tdc_timing_offsets_psc.txt
ccdb add /PHOTON_BEAM/pair_spectrometer/coarse/adc_timing_offsets -v default -r ${RunNo}-${RunNo} ${Dir}/adc_timing_offsets_psc.txt
ccdb add /PHOTON_BEAM/pair_spectrometer/fine/adc_timing_offsets -v default -r ${RunNo}-${RunNo} ${Dir}/adc_timing_offsets_ps.txt
