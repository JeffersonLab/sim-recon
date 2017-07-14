void call_bcal_hadronic_eff(TString file, TString option)
{
// issue the tree->Process, so that it can be run from the command line
//
    cout << "call_bcal_hadronic_eff, file=" << file << " option=" << option << endl;
    bcal_hadronic_eff->Process(file,option);
}

