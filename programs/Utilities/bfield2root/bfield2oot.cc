// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "TRACKING/DMagneticFieldMap.h"
#include "JANA/JParameterManager.h"

#include <TFile.h>
#include <TTree.h>

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	float val[10];

	// open ROOT file
	TFile* ROOTfile = new TFile("bfield.root","RECREATE","Produced by bfield2root");
	ROOTfile->SetCompressionLevel(6);
	cout<<"Opened ROOT file \"bfield.root\""<<endl;

	// Create Magnetic Field Map
	new JParameterManager(); // normally, this is created by JApplication!
	DMagneticFieldMap *bfield = new DMagneticFieldMap();

	// Create Tree
	TTree *tree = new TTree("Bfield","Magnetic Field");
	tree->Branch("B",val,"x/F:y:z:r:phi:Bx:By:Bz:Br:Bphi");
	
	// Get the field map
	const DBfieldPoint_t* Bmap;
	int Npoints;
	bfield->GetTable(Bmap, Npoints);
	
	double BMAP_Z_OFFSET;
	gPARMS->GetParameter("GEOM:BMAP_Z_OFFSET", BMAP_Z_OFFSET);
	

	// Loop over all points in the map
	const DBfieldPoint_t* B = Bmap;
	for(int i=0; i<Npoints; i++, B++){
		double x = B->x*2.54;
		double y = B->y*2.54;
		double z = (B->z-BMAP_Z_OFFSET)*2.54;
		double Bx = B->Bx;
		double By = B->By;
		double Bz = B->Bz;
		double r = sqrt(x*x+y*y);
		double phi = atan2(y,x);
		double Br = Bx;
		double Bphi = By;

		val[0] = x;
		val[1] = y;
		val[2] = z;
		val[3] = r;
		val[4] = phi;
		val[5] = Bx;
		val[6] = By;
		val[7] = Bz;
		val[8] = Br;
		val[9] = Bphi;
		
		tree->Fill();
	}
		
	ROOTfile->Write();
	delete ROOTfile;
	cout<<endl<<"Closed ROOT file"<<endl;
		
	return 0;
}

