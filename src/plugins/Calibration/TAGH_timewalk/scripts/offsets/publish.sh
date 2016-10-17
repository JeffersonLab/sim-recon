RunNo=$1
RunNo_min=$2
RunNo_max=$3
Dir=OUTPUT/${RunNo}/offsets
ccdb add /PHOTON_BEAM/hodoscope/base_time_offset -v default -r ${RunNo_min}-${RunNo_max} ${Dir}/base_time_offset.txt
ccdb add /PHOTON_BEAM/hodoscope/tdc_time_offsets -v default -r ${RunNo_min}-${RunNo_max} ${Dir}/tdc_time_offsets.txt
ccdb add /PHOTON_BEAM/hodoscope/fadc_time_offsets -v default -r ${RunNo_min}-${RunNo_max} ${Dir}/fadc_time_offsets.txt
