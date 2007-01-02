// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

#include "TRACKING/DMagneticFieldMapGlueX.h"
#include "JANA/JParameterManager.h"

#include <TFile.h>
#include <TTree.h>

extern "C"{
 void gufld_(float*, float*); // FORTRAN
};

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
	DMagneticFieldMap *bfield = new DMagneticFieldMapGlueX();

	// Create Trees
	TTree *gtree = new TTree("gfield","HDGeant Magnetic Field");
	gtree->Branch("B",val,"x/F:y:z:r:phi:Bx:By:Bz:Br:Bphi");

	TTree *btree = new TTree("bfield","DANA Magnetic Field");
	btree->Branch("B",val,"x/F:y:z:r:phi:Bx:By:Bz:Br:Bphi");
	
	// Space grid
	for(float r=0.0; r<65.0; r+=2.0){
		for(float phi=0.0; phi<2.0*M_PI; phi+=M_PI/50.0){
			for(float z=-50.0; z<650.0; z+=5.0){
				float u[3], B[3];
				u[0] = r*cos(phi);
				u[1] = r*sin(phi);
				u[2] = z;
				gufld_(u,B);
				//B[0]=B[1]=B[2] = 0.0;
				
				// Geant uses kGauss. Convert to Tesla
				B[0]/=10.0;
				B[1]/=10.0;
				B[2]/=10.0;
			
				double Br = (B[0]*u[0] + B[1]*u[1])/(r+1.0E-20);
				double Bphi = (B[1]*u[0] - B[0]*u[1])/(r+1.0E-20);

				val[0] = u[0];
				val[1] = u[1];
				val[2] = u[2];
				val[3] = r;
				val[4] = phi;
				val[5] = B[0];
				val[6] = B[1];
				val[7] = B[2];
				val[8] = Br;
				val[9] = Bphi;

				gtree->Fill();
				
				// Access DANA field map at the same point
				double x,y,z,Bx,By,Bz;
				bfield->GetField(x=val[0], y=val[1], z=val[2], Bx, By, Bz);
				B[0] = Bx;
				B[1] = By;
				B[2] = Bz;

				Br = (B[0]*u[0] + B[1]*u[1])/(r+1.0E-20);
				Bphi = (B[1]*u[0] - B[0]*u[1])/(r+1.0E-20);

				val[0] = u[0];
				val[1] = u[1];
				val[2] = u[2];
				val[3] = r;
				val[4] = phi;
				val[5] = B[0];
				val[6] = B[1];
				val[7] = B[2];
				val[8] = Br;
				val[9] = Bphi;

				btree->Fill();
			}
		}
	}
		
	ROOTfile->Write();
	delete ROOTfile;
	cout<<endl<<"Closed ROOT file"<<endl;
		
	return 0;
}

