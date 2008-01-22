// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "HDGEOMETRY/DMagneticFieldMapGlueX.h"
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
	DMagneticFieldMapGlueX *bfield = new DMagneticFieldMapGlueX();

	// Create Tree
	TTree *tree = new TTree("Bfield","Magnetic Field");
	tree->Branch("B",val,"x/F:y:z:r:phi:Bx:By:Bz:Br:Bphi");
	
	// Get the field map
	const DMagneticFieldMapGlueX::DBfieldPoint_t* Bmap;
	int Npoints;
	bfield->GetTable(Bmap, Npoints);

#if 0
	// Space grid
	for(float r=0.0; r<100.0; r+=2.0){
		for(float phi=0.0; phi<2.0*M_PI; phi+=M_PI/50.0){
			for(float z=-1000.0; z<1000.0; z+=5.0){
				float u[3], B[3];
				u[0] = r*cos(phi);
				u[1] = r*sin(phi);
				u[2] = z;
				gufld_(u,B);
				//B[0]=B[1]=B[2] = 0.0;
			
				double Br = (B[0]*u[0] + B[1]*u[1])/(r+1.0E-20);
				double Bphi = (B[1]*u[0] - B[0]*u[1])/(r+1.0E-20);
				double Btot = sqrt(Br*Br + B[2]*B[2]);
				if(Btot==0.0)continue;

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

				tree->Fill();
			}
		}
	}

#else

	double BMAP_Z_OFFSET;
	gPARMS->GetParameter("GEOM:BMAP_Z_OFFSET", BMAP_Z_OFFSET);
	

	// Loop over all points in the map
	const DMagneticFieldMapGlueX::DBfieldPoint_t* B = Bmap;
	for(int i=0; i<Npoints; i++, B++){
		double x = B->x;
		double y = B->y;
		double z = B->z;
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
#endif
		
	ROOTfile->Write();
	delete ROOTfile;
	cout<<endl<<"Closed ROOT file"<<endl;
		
	return 0;
}

