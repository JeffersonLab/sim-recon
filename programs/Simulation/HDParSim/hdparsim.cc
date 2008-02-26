
// This is s simple example program that shows how to use a DTrackResolution
// class to apply detector derived resolutions to charged tracks.
// There are detailed comments that should step you through the program
// with enough information to understand it and allow you to modify it to
// your own needs.


#include <iostream>
using namespace std;

#include <particleType.h>

#include "DTrackingResolutionGEANT.h"

//------------
// main
//------------
int main(int narg, char *argv[])
{
	// Create a DTrackingResolution object. Specifically,
	// we create a DTrackingResolutionGEANT object so as
	// to use the resolutions derived from the GEANT-based
	// simulations.
	DTrackingResolution *res = new DTrackingResolutionGEANT();

	// Open a root file to hold our output histogram(s)
	TFile *f = new TFile("hdparsim.root","RECREATE");
	if(!f->IsOpen()){
		cerr<<"Unable to open ROOT output file!"<<endl;
		return -1;
	}
	
	// Create a histogram of dp over p vs. p
	TH2D *dp_over_p_vs_p = new TH2D("dp_over_p_vs_p","#deltap/p vs. p", 100, 0.0, 7.0, 100, -0.2, 0.2);
	
	// Randomly pick 1E6, points in p and smear the momentum
	// so we can calculate dp/p. Randomly pick phi also, but
	// throw everything at theta=30 degrees.
	TRandom3 rnd;
	double theta = 30.0/57.3;
	for(int i=0; i<1E6; i++){
		
		// Randomly pick momentum and phi angle.
		// These are our "true" thrown values.
		double ptot = 0.2 + rnd.Rndm()*6.8; // pick a momentum between 0.2 and 7
		double phi = rnd.Rndm()*2.0*M_PI; // pick a phi between 0 and 2pi
		TVector3 mom;
		mom.SetMagThetaPhi(ptot, theta, phi);
		
		// Smear momentum. On input, "mom" has the true values and on
		// output, it contains the smeared values.
		res->Smear(PiPlus, mom);
		
		// Fill histogram
		dp_over_p_vs_p->Fill(ptot, (ptot-mom.Mag())/mom.Mag());
	}

	// Close ROOT output file
	f->Write();
	delete f;

	return 0;
}

