RunNo=$1
OutputDir=OUTPUT/$RunNo
InputFile=$2
cp -p ${InputFile} hd_root.root
root -b -q 'gaussian_fits.C("hd_root.root",true)'
root -b -q 'timewalk_fits.C("gaussian-fits-csv")'
rm -f hd_root.root
mkdir -p $OutputDir
mv *.txt $OutputDir; mv fits_* $OutputDir; mv parms_timewalk $OutputDir
mv gaussian-fits-csv $OutputDir; mv overall_gaussian_fit.gif $OutputDir
