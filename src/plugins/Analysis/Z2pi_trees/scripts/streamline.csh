#!/bin/csh -f
set echo
#
# streamline.csh
# Elton Smith. Sep 18, 2017. 
# Streamline instructions to process the output of mcsmear output files through amplitude analysis.
#
source ~/MC_environment.csh
cd /work/halld/home/elton/gen_2pi_primakoff_signal/Z2pi_trees
unset run
unset base

if( $#argv == 0 ) then
    set run = "031000"
    set base = "signal"
else if ( $#argv == 1) then
    set run = $1
    set base = "signal"
else if ( $#argv == 1) then
    set run = $1
    set base = $2
else
    echo "Too many arguments=" $1 $2 $3
endif

echo " run =" $run
echo " base=" $base

# gen_2pi_primakoff> gen_2pi_primakoff -c gen_2pi_primakoff_${base}.cfg -o tree_gen_2pi_primakoff_${base}.root -hd gen_2pi_primakoff_${base}.hddm -a 5.5 -b 6.0 -p 6.0 -m 11.6 -n 100000 -r 31000

# Here are instructions for processing MC smeared output files / or data

# hd_root -PPLUGINS=monitoring_hists,Z2pi_trees -PNTHREADS=4 ../../gen_2pi_primakoff_${base}/hddm/dana_rest_gen_2pi_primakoff_${run}_*.hddm -o hd_root_Z2pi_trees_${base}.root
# mv tree_Z2pi_trees.root tree_hd_root_Z2pi_trees_${base}.root
root -b -q tree_hd_root_Z2pi_trees_${base}.root 'call_DSelector.C("DSelector_Z2pi_trees.C+")'
mv DSelector_Z2pi_trees.root DSelector_Z2pi_trees_${base}.root
mv tree_DSelector_Z2pi_trees.root tree_DSelector_Z2pi_trees_${base}.root
root -l -q  plot_Z2pi_trees.C\(\"DSelector_Z2pi_trees_${base}\"\)
tree_to_amptools tree_DSelector_Z2pi_trees_${base}.root Z2pi_trees_Tree
mv AmpToolsInputTree.root tree_DSelector_Z2pi_trees_${base}_amptools.root

# Now repeat for flat distribution
set base = "flat"

# gen_2pi_primakoff> gen_2pi_primakoff -c gen_2pi_primakoff_${base}.cfg -o tree_gen_2pi_primakoff_${base}.root -hd gen_2pi_primakoff_${base}.hddm -a 5.5 -b 6.0 -p 6.0 -m 11.6 -n 1000000 -r 31000

# hd_root -PPLUGINS=monitoring_hists,Z2pi_trees -PNTHREADS=4 ../../gen_2pi_primakoff_${base}/hddm/dana_rest_gen_2pi_primakoff_${run}_*.hddm -o hd_root_Z2pi_trees_${base}.root
# mv tree_Z2pi_trees.root tree_hd_root_Z2pi_trees_${base}.root
root -b -q tree_hd_root_Z2pi_trees_${base}.root 'call_DSelector.C("DSelector_Z2pi_trees.C+")'
mv DSelector_Z2pi_trees.root DSelector_Z2pi_trees_${base}.root
mv tree_DSelector_Z2pi_trees.root tree_DSelector_Z2pi_trees_${base}.root
root -l -q  plot_Z2pi_trees.C\(\"DSelector_Z2pi_trees_${base}\"\)
tree_to_amptools tree_DSelector_Z2pi_trees_${base}.root Z2pi_trees_Tree
mv AmpToolsInputTree.root tree_DSelector_Z2pi_trees_${base}_amptools.root


fit -c fit_2pi_primakoff.cfg
cp twopi_primakoff.fit twopi_primakoff_DSelect.fit
twopi_plotter_primakoff twopi_primakoff_DSelect.fit -o twopi_primakoff_DSelect.root
mv twopi_fitPars.txt twopi_primakoff_DSelect.fit2
# root -q -l twopi_primakoff.C

unset echo

