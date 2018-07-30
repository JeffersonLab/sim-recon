void call_DSelector (TString file)
{
// issue the tree->Process, so that it can be run from the command line
//
cout << "call_DSelector: file=" << file << endl;
gROOT->LoadMacro("$ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
Z2pi_trees_Tree->Process("DSelector_Z2pi_trees.C+");
}
