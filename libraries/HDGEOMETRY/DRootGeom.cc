// $Id$
//
//    File: DRootGeom.h
// Created: Fri Feb 13 08:43:39 EST 2009
// Creator: zihlmann
//

#include "DRootGeom.h"
#include "hddsroot.h"

using namespace std;


//---------------------------------
// DRootGeom    (Constructor)
//---------------------------------
DRootGeom::DRootGeom()
{
	pthread_mutexattr_init(&mutex_attr);
	pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_RECURSIVE);
	pthread_mutex_init(&mutex, &mutex_attr);

	DRGeom = hddsroot();

	// Fill average materials table
	cout<<"Filling material table ..."; cout.flush();
	Nr = 125;
	Nz = 750;
	double rmin = 0.0;
	double rmax = 125.0;
	double zmin = -100.0;
	double zmax = 650.0;
	r0 = rmin;
	z0 = zmin;
	dr = (rmax-rmin)/(double)Nr;
	dz = (zmax-zmin)/(double)Nz;
	MatTable = new VolMat*[Nr];
	buff = new VolMat[Nr*Nz];
	for(int ir=0; ir<Nr; ir++){
		double r = r0 + (double)ir*dr;
		MatTable[ir]=&buff[ir*Nz];
		for(int iz=0; iz<Nz; iz++){
			double z = z0 + (double)iz*dz;
			
			// Loop over points in phi, r, and z and add up material
			int n_r = 2;
			int n_z = 2;
			int n_phi = 20;
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
						FindMat(pos, density, A, Z, radlen);
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
}

//---------------------------------
// DRootGeom    (Destructor)
//---------------------------------
DRootGeom::~DRootGeom()
{
  delete DRGeom;
}

TGeoNode* DRootGeom::FindNode(double *x)
{
	pthread_mutex_lock(&mutex);

  TGeoNode *cnode = DRGeom->FindNode(x[0],x[1],x[2]);
  
	pthread_mutex_unlock(&mutex);

  return cnode;  

}

//---------------------------------
// FindVolume
//---------------------------------
TGeoVolume* DRootGeom::FindVolume(double *x)
{

	pthread_mutex_lock(&mutex);
  TGeoNode *cnode = DRGeom->FindNode(x[0],x[1],x[2]);
  TGeoVolume *cvol = cnode->GetVolume();
	pthread_mutex_unlock(&mutex);

  return cvol;  

}

//---------------------------------
// FindMat
//---------------------------------
jerror_t DRootGeom::FindMat(DVector3 pos,double &density, double &A, double &Z, double &RadLen) const
{
	return FindMatLL(pos, density, A, Z, RadLen);
}

//---------------------------------
// FindMatTable
//---------------------------------
jerror_t DRootGeom::FindMatTable(DVector3 pos,double &density, double &A, double &Z, double &RadLen) const
{
	// For now, this just finds the bin in the material map the given position is in
	// (i.e. no interpolation )
	double r = pos.Perp();
	double z = pos.Z();
	int ir = floor((r-r0)/dr);
	int iz = floor((z-z0)/dz);
	if(ir<0 || ir>=Nr || iz<0 || iz>=Nz){
		A = Z = density = RadLen = 0.0;
		return RESOURCE_UNAVAILABLE;
	}
	
	const VolMat &mat = MatTable[ir][iz];
	A = mat.A;
	Z = mat.Z;
	density = mat.Density;
	RadLen = mat.RadLen;

	return NOERROR;
}

//---------------------------------
// FindMatLL
//---------------------------------
jerror_t DRootGeom::FindMatLL(DVector3 pos,double &density, double &A, double &Z,
			double &RadLen) const{
  density=RadLen=A=Z=0.;

	pthread_mutex_lock(const_cast<pthread_mutex_t*>(&mutex));

  TGeoNode *cnode = DRGeom->FindNode(pos.X(),pos.Y(),pos.Z());
  if (cnode==NULL){
    _DBG_<<"Missing cnode at position (" <<pos.X()<<","<<pos.Y()<<","
	 <<pos.Z()<<")"<<endl;
    return RESOURCE_UNAVAILABLE;
  }
  TGeoVolume *cvol = cnode->GetVolume();
  if (cvol==NULL){
    _DBG_<<"Missing cvol" <<endl;
    return RESOURCE_UNAVAILABLE;
  }
  TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial();

  density=cmat->GetDensity();
  RadLen=cmat->GetRadLen();
  A=cmat->GetA();
  Z=cmat->GetZ();
  if (A<1.){ // try to prevent division by zero problems
    A=1.;
  } 

	pthread_mutex_unlock(const_cast<pthread_mutex_t*>(&mutex));

  return NOERROR;
}


//---------------------------------
// FindMat
//---------------------------------
struct VolMat DRootGeom::FindMat(double *x)
{

	pthread_mutex_lock(&mutex);
  TGeoNode *cnode = DRGeom->FindNode(x[0],x[1],x[2]);
  TGeoVolume *cvol = cnode->GetVolume();
  TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial();

  Current_Node = cnode;
  Current_Volume = cvol;
  Current_Material = cmat;
  Mat_Index = cmat->GetIndex();

  Mat.A = cmat->GetA();
  Mat.Z = cmat->GetZ();
  Mat.Density = cmat->GetDensity();
  Mat.RadLen = cmat->GetRadLen();

	pthread_mutex_unlock(&mutex);

  //cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
  //cout<<DRGeom->FindNode(x[0],x[1],x[2])->GetVolume()->GetName()<<endl; 

  return Mat;  

}


