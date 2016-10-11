RunNo=$1
OutputDir=OUTPUT/${RunNo}
InputFile=$2
cp -p ${InputFile} hd_root.root
mkdir -p offsets
ccdb dump /PHOTON_BEAM/hodoscope/base_time_offset:${RunNo} > offsets/base_time_offset_ccdb.txt
ccdb dump /PHOTON_BEAM/hodoscope/tdc_time_offsets:${RunNo} > offsets/tdc_time_offsets_ccdb.txt
ccdb dump /PHOTON_BEAM/hodoscope/fadc_time_offsets:${RunNo} > offsets/fadc_time_offsets_ccdb.txt
root -b -q 'fits.C("hd_root.root",true)'
root -b -q 'offsets.C("fits-csv")'
rm -f hd_root.root
mkdir -p $OutputDir
mv offsets $OutputDir; mv fits-csv $OutputDir; mv fits_* $OutputDir
