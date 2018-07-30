#!/bin/csh -f
set echo
#
# streamline.csh
# Elton Smith. Sep 18, 2017. 
# Streamline instructions to process the output of mcsmear output files through amplitude analysis.
#
source ~/MC_environment.csh
unset run
unset base
unset maxev

if( $#argv == 0 ) then
    set run = "031000"
    set base = "NoBeamPipe"
    set maxev = 100000
else if ( $#argv == 1) then
    set run = $1
    set base = "NoBeamPipe"
    set maxev = 100000
else if ( $#argv == 2) then
    set run = $1
    set base = $2
    set maxev = 100000
else if ( $#argv == 3) then
    set run = $1
    set base = $2
    set maxev = $3
else
    echo "Too many arguments=" $1 $2 $3 $4
endif

cd /work/halld/home/elton/gen_2pi_primakoff_${base}/Z2pi_trees

echo " run =" $run
echo " base=" $base
echo " maxev=" $maxev

# Here are instructions for processing MC smeared output files / or data

hd_root -PPLUGINS=monitoring_hists,Z2pi_trees -PNTHREADS=4  -PEVENTS_TO_KEEP=${maxev} ../../gen_2pi_primakoff_${base}/hddm/dana_rest_gen_2pi_primakoff_${base}_signal_${run}_*.hddm -o hd_root_Z2pi_trees_${base}_signal_${maxev}.root
mv tree_Z2pi_trees.root tree_hd_root_Z2pi_trees_${base}_signal_${maxev}.root
root -b -q tree_hd_root_Z2pi_trees_${base}_signal_${maxev}.root 'call_DSelector.C("DSelector_Z2pi_trees.C+")'
mv DSelector_Z2pi_trees.root DSelector_Z2pi_trees_${base}_signal_${maxev}.root
mv tree_DSelector_Z2pi_trees.root tree_DSelector_Z2pi_trees_${base}_signal_${maxev}.root
root -l -q  plot_Z2pi_trees.C\(\"DSelector_Z2pi_trees_${base}_signal_${maxev}\"\)
tree_to_amptools tree_DSelector_Z2pi_trees_${base}_signal_${maxev}.root Z2pi_trees_Tree
mv AmpToolsInputTree.root tree_DSelector_Z2pi_trees_${base}_signal_${maxev}_amptools.root

# Now repeat for flat distribution. Also need the generated flat distributions.
# set base = "NoBeamPipe"

gen_2pi_primakoff -c gen_2pi_primakoff_flat.cfg -o tree_gen_2pi_primakoff_flat_${maxev}.root -hd gen_2pi_primakoff_flat.hddm -a 5.5 -b 6.0 -p 6.0 -m 11.6 -n ${maxev} -r ${run}

hd_root -PPLUGINS=monitoring_hists,Z2pi_trees -PNTHREADS=4 -PEVENTS_TO_KEEP=${maxev} ../../gen_2pi_primakoff_${base}/hddm/dana_rest_gen_2pi_primakoff_${base}_flat_${run}_*.hddm -o hd_root_Z2pi_trees_${base}_flat_${maxev}.root
mv tree_Z2pi_trees.root tree_hd_root_Z2pi_trees_${base}_flat_${maxev}.root
root -b -q tree_hd_root_Z2pi_trees_${base}_flat_${maxev}.root 'call_DSelector.C("DSelector_Z2pi_trees.C+")'
mv DSelector_Z2pi_trees.root DSelector_Z2pi_trees_${base}_flat_${maxev}.root
mv tree_DSelector_Z2pi_trees.root tree_DSelector_Z2pi_trees_${base}_flat_${maxev}.root
root -l -q  plot_Z2pi_trees.C\(\"DSelector_Z2pi_trees_${base}_flat_${maxev}\"\)
tree_to_amptools tree_DSelector_Z2pi_trees_${base}_flat_${maxev}.root Z2pi_trees_Tree
mv AmpToolsInputTree.root tree_DSelector_Z2pi_trees_${base}_flat_${maxev}_amptools.root


fit -c fit_2pi_primakoff_${maxev}.cfg
cp twopi_primakoff.fit twopi_primakoff_DSelect_${base}_${maxev}.fit
twopi_plotter_primakoff twopi_primakoff_DSelect_${base}_${maxev}.fit -o twopi_primakoff_DSelect_${base}_${maxev}.root
mv twopi_fitPars.txt twopi_primakoff_DSelect_${base}_${maxev}.fit2
root -q -l twopi_primakoff.C\(\"twopi_primakoff_DSelect_${base}_${maxev}\",${maxev}\)

unset echo

