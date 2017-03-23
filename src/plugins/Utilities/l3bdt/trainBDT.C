//
// May 12, 2016  D. Lawrence
//
// This macro is used to train a BDT to be used in
// the Level-3 trigger. It assumes there already
// exists a file l3bdt.root that was created using
// the l3bdt plugin like this:
//
//  hd_root -o l3bdt.root -PPLUGINS=L3BDTtree file.evio
//

#include <iostream>

#include<TMVA/Factory.h>

void trainBDT(void)
{
	// Open input file and get tree
	TFile *infile = new TFile("l3bdt.root");
	TTree *l3tree = (TTree*)infile->Get("l3tree");
	if(l3tree == NULL){
		cout << "Couldn't open \"l3bdt.root\"!" << endl;
		return;
	}

	// Open output root file (for TMVA)
	TFile *outfile = new TFile("l3BDT_out.root", "RECREATE");
	TMVA::Factory *fac = new TMVA::Factory("L3",outfile,"");

	// Specify input tree that contains both signal and background 
	TCut signalCut("Evisible>=4.0 && ((Npshits+Npschits)<=1)");
	TCut backgroundCut("Evisible<4.0 && ((Npshits+Npschits)<=1)");
	fac->SetInputTrees(l3tree, signalCut, backgroundCut);

	// Add variables
	fac->AddVariable("Nstart_counter",      'I');
	fac->AddVariable("Ntof",                'I');
	fac->AddVariable("Nbcal_points",        'I');
	fac->AddVariable("Nbcal_clusters",      'I');
	fac->AddVariable("Ebcal_points",        'F');
	fac->AddVariable("Ebcal_clusters",      'F');
	fac->AddVariable("Nfcal_clusters",      'I');
	fac->AddVariable("Efcal_clusters",      'F');
	fac->AddVariable("Ntrack_candidates",   'I');
	fac->AddVariable("Ptot_candidates",     'F');

	TCut preSelectCut("");
	fac->PrepareTrainingAndTestTree(preSelectCut,"");
	fac->BookMethod(TMVA::Types::kBDT, "BDT", "");

	fac->TrainAllMethods();
	fac->TestAllMethods();
	fac->EvaluateAllMethods();

	delete fac;

	outfile->Close();
	delete outfile;
}

