// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DMagneticFieldMap.h"

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
	cout<<"Opened ROOT file \"bfield.root\""<<endl;

	// Create Magnetic Field Map
	DMagneticFieldMap *bfield = new DMagneticFieldMap(41, 251);
	bfield->readMap();

	// Create Tree
	TTree *B = new TTree("Bfield","Magnetic Field");
	B->Branch("B",val,"x/F:y:z:r:phi:Bx:By:Bz:Br:Bphi");
	
	// Get Dimensions of field map
	int rDim = bfield->Get_rDim();
	int zDim = bfield->Get_zDim();
	int phiDim = bfield->Get_phiDim();
	double rMin = bfield->Get_rMin();
	double rMax = bfield->Get_rMax();
	double zMin = bfield->Get_zMin();
	double zMax = bfield->Get_zMax();
	const double *map = bfield->Get_map();
	
	double rstep = (rMax-rMin)/(double)(rDim-1);
	double zstep = (zMax-zMin)/(double)(zDim-1);
	double phistep = 2.0*M_PI/(double)(phiDim-1);
	
	float r = rMin;
	float z = zMin;
	float phi = 0.0;
	for(int ir=0;ir<rDim;ir++, r+=rstep){
		z=zMin;
		for(int iz=0;iz<zDim;iz++, z+=zstep){
			phi=0.0;
			for(int iphi=0;iphi<=phiDim;iphi++, phi+=phistep){
				float x = r*cos(phi);
				float y = r*sin(phi);

				int ind = bfield->serialize(ir,iz);
				float Br = map[ind+0];
				float Bphi = map[ind+1];
				float Bz = map[ind+2];
				
				val[0] = x;
				val[1] = y;
				val[2] = z;
				val[3] = r;
				val[4] = phi;
				val[5] = Br*cos(phi);  // needs to add Bphi component!!
				val[6] = Br*sin(phi);  // needs to add Bphi component!!
				val[7] = Bz;
				val[8] = Br;
				val[9] = Bphi;
				
				cout<<"Filling B (Bz="<<Bz<<")"<<endl;
				B->Fill();
			}
		}
	}
	
	ROOTfile->Write();
	delete ROOTfile;
	cout<<endl<<"Closed ROOT file"<<endl;
		
	return 0;
}

