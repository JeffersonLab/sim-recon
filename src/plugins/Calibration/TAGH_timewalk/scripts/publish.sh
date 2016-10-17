RunNo=$1
RunNo_min=$2
RunNo_max=$3
Dir=OUTPUT/${RunNo}
ccdb add /PHOTON_BEAM/hodoscope/tdc_timewalk -v default -r ${RunNo_min}-${RunNo_max} ${Dir}/tdc_timewalk.txt
