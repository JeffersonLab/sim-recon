// Author: David Lawrence  June 19, 2009
//
//
// mkMaterialMap.cc
//
#include <fstream>

#include "DANA/DApplication.h"
using namespace std;

#include <DVector3.h>
#include <HDGEOMETRY/DRootGeom.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH2D.h>

void Usage(JApplication &app);


//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	DApplication *dapp = new DApplication(narg, argv);

	DRootGeom *g = new DRootGeom(dapp);

	// Fill average materials table
	cout<<"Filling material table ..."; cout.flush();
	int Nr = 125;
	int Nz = 750;
	double rmin = 0.0;
	double rmax = 125.0;
	double zmin = -100.0;
	double zmax = 650.0;
	double r0 = rmin;
	double z0 = zmin;
	double dr = (rmax-rmin)/(double)Nr;
	double dz = (zmax-zmin)/(double)Nz;
	VolMat **MatTable = new VolMat*[Nr];
	VolMat *buff = new VolMat[Nr*Nz];
	for(int ir=0; ir<Nr; ir++){
		double r = r0 + (double)ir*dr;
		MatTable[ir]=&buff[ir*Nz];
		for(int iz=0; iz<Nz; iz++){
			double z = z0 + (double)iz*dz;
			
			// Loop over points in phi, r, and z and add up material
			int n_r = 5;
			int n_z = 5;
			int n_phi = 60;
			double d_r = dr/(double)n_r;
			double d_z = dz/(double)n_z;
			double d_phi = 2.0*M_PI/(double)n_phi;
			VolMat avg_mat={0.0, 0.0, 0.0, 0.0};
			for(int i_r=0; i_r<n_r; i_r++){
				double my_r = r - dr/2.0 + (double)i_r*d_r;
				for(int i_z=0; i_z<n_z; i_z++){
					double my_z = z - dz/2.0 + (double)i_z*d_z;
					for(int i_phi=0; i_phi<n_phi; i_phi++){
						double my_phi = (double)i_phi*d_phi;

						DVector3 pos(my_r*cos(my_phi), my_r*sin(my_phi), my_z);
						double A, Z, density, radlen;
						g->FindMatLL(pos, density, A, Z, radlen);
						avg_mat.A += A;
						avg_mat.Z += Z;
						avg_mat.Density += density;
						avg_mat.RadLen += radlen;
					}
				}
			}

			// Divide by number of points to get averages
			avg_mat.A /= (double)(n_r*n_z*n_phi);
			avg_mat.Z /= (double)(n_r*n_z*n_phi);
			avg_mat.Density /= (double)(n_r*n_z*n_phi);
			avg_mat.RadLen /= (double)(n_r*n_z*n_phi);
			
			MatTable[ir][iz] = avg_mat;
		}
		cout<<"\r Filling Material table ... "<<100.0*(double)ir/(double)Nr<<"%       ";cout.flush();
	}
	cout <<"Done"<<endl;
	
	// Write table to file
	cout<<"Writing material table to file ..."<<endl;
	ofstream of("material_map"); 
	of<<"#  r	z	A	Z	density radlen"<<endl;
	for(int ir=0; ir<Nr; ir++){
		double r = r0 + (double)ir*dr;
		for(int iz=0; iz<Nz; iz++){
			double z = z0 + (double)iz*dz;
			
			VolMat &mat = MatTable[ir][iz];
			
			of<<r<<"\t"<<z<<"\t"<<mat.A<<"\t"<<mat.Z<<"\t"<<mat.Density<<"\t"<<mat.RadLen<<endl;
		}
	}
	of.close();
	
	return 0;
}

//-----------
// Usage
//-----------
void Usage(JApplication &app)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"    mkMaterialMap [options] source1 source2 source3 ..."<<endl;
	cout<<endl;
	app.Usage();
	cout<<endl;
	
	exit(0);
}

