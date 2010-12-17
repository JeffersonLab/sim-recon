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
DRootGeom::DRootGeom(JApplication *japp)
{
	pthread_mutexattr_init(&mutex_attr);
	pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_RECURSIVE);
	pthread_mutex_init(&mutex, &mutex_attr);

	table_initialized = false;
	DRGeom = NULL;

	int runnumber = 1;
	jcalib = japp->GetJCalibration(runnumber);
	if(!jcalib){
		_DBG_<<"Unable to get JCalibration object!"<<endl;
		exit(-1);
	}

	JParameterManager *jparms = japp->GetJParameterManager();
	if(!jparms){
		_DBG_<<"Unable to get JParameterManager object!"<<endl;
		exit(-1);
	}

}

//---------------------------------
// DRootGeom    (Destructor)
//---------------------------------
DRootGeom::~DRootGeom()
{
  delete DRGeom;
}

//---------------------------------
// ReadMap
//---------------------------------
int DRootGeom::ReadMap(string namepath, int runnumber)
{
	/// Read the magnetic field map in from the calibration database.
	/// This will read in the map and figure out the number of grid
	/// points in each direction (x,y, and z) and the range in each.
	/// The gradiant of the field is calculated for all but the most
	/// exterior points and saved to use in later calls to GetField(...).

	// Read in map from calibration database. This should really be
	// built into a run-dependent geometry framework, but for now
	// we do it this way. 
	if(!jcalib)return 0;
	
	cout<<"Reading Material map from "<<namepath<<" ..."<<endl;
	vector< vector<float> > Mmap;
	jcalib->Get(namepath, Mmap);
	cout<<Mmap.size()<<" entries found (";
	if(Mmap.size()<1){
		cout<<")"<<endl;
		return Mmap.size();
	}
	
	// The map should be on a grid with equal spacing in r, and z.
	// Here we want to determine the number of points in each of these
	// dimensions and the range. 
	// The easiest way to do this is to use a map<float, int> to make a
	// histogram of the entries by using the key to hold the extent
	// so that the number of entries will be equal to the number of
	// different values.
	map<float, int> rvals;
	map<float, int> zvals;
	double rmin, zmin, rmax, zmax;
	rmin = zmin = 1.0E6;
	rmax = zmax = -1.0E6;
	for(unsigned int i=0; i<Mmap.size(); i++){
		vector<float> &a = Mmap[i];
		float &r = a[0];
		float &z = a[1];
		
		rvals[r] = 1;
		zvals[z] = 1;
		if(r<rmin)rmin=r;
		if(z<zmin)zmin=z;
		if(r>rmax)rmax=r;
		if(z>zmax)zmax=z;
	}
	Nr = rvals.size();
	Nz = zvals.size();
	r0 = rmin;
	z0 = zmin;
	dr = (rmax-rmin)/(double)(Nr-1);
	dz = (zmax-zmin)/(double)(Nz-1);
	cout<<" Nr="<<Nr;
	cout<<" Nz="<<Nz;
	cout<<" )  at 0x"<<hex<<(unsigned long)this<<dec<<endl;
	
	// Create 2D array to hold table in class
	MatTable = new VolMat*[Nr];
	buff = new VolMat[Nr*Nz];
	for(int ir=0; ir<Nr; ir++){
		MatTable[ir]=&buff[ir*Nz];
	}
	
	// Fill table
	for(unsigned int i=0; i<Mmap.size(); i++){
		vector<float> &a = Mmap[i];
		float &r = a[0];
		float &z = a[1];
		int ir = (int)floor((r - r0)/dr);
		int iz = (int)floor((z - z0)/dz);
		if(ir<0 || ir>=Nr){_DBG_<<"ir out of range: ir="<<ir<<"  Nr="<<Nr<<endl; continue;}
		if(iz<0 || iz>=Nz){_DBG_<<"iz out of range: iz="<<iz<<"  Nz="<<Nz<<endl; continue;}
		VolMat &mat = MatTable[ir][iz];
		mat.A = a[2];
		mat.Z = a[3];
		mat.Density = a[4];
		mat.RadLen = a[5];
		mat.rhoZ_overA = a[6];
		mat.rhoZ_overA_logI = a[7];
	}
	
	
	return Mmap.size();
}

//---------------------------------
// InitTable
//---------------------------------
void DRootGeom::InitTable(void)
{
	pthread_mutex_lock(&mutex);
	
	if(table_initialized){
		// Don't initialize table twice!
		pthread_mutex_unlock(&mutex);
		return;
	}
	
	// Fill average materials table
	cout<<"Filling material table"; cout.flush();
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
			int n_r = 3;
			int n_z = 2;
			int n_phi = 20;
			double d_r = dr/(double)n_r;
			double d_z = dz/(double)n_z;
			double d_phi = 2.0*M_PI/(double)n_phi;
			VolMat avg_mat={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

			for(int i_r=0; i_r<n_r; i_r++){
				double my_r = r - dr/2.0 + (double)i_r*d_r;
				for(int i_z=0; i_z<n_z; i_z++){
					double my_z = z - dz/2.0 + (double)i_z*d_z;
					for(int i_phi=0; i_phi<n_phi; i_phi++){
						double my_phi = (double)i_phi*d_phi;

						DVector3 pos(my_r*cos(my_phi), my_r*sin(my_phi), my_z);
						double A, Z, density, radlen,LnI;

						FindMatLL(pos, density, A, Z, radlen,LnI);

						double rhoZ_overA = density*Z/A;
						double rhoZ_overA_logI = rhoZ_overA*LnI;

						avg_mat.A += A;
						avg_mat.Z += Z;
						avg_mat.Density += density;
						avg_mat.RadLen += radlen;
						avg_mat.rhoZ_overA += rhoZ_overA;
						avg_mat.rhoZ_overA_logI += rhoZ_overA_logI;
					}
				}
			}

			// Divide by number of points to get averages
			avg_mat.A /= (double)(n_r*n_z*n_phi);
			avg_mat.Z /= (double)(n_r*n_z*n_phi);
			avg_mat.Density /= (double)(n_r*n_z*n_phi);
			avg_mat.RadLen /= (double)(n_r*n_z*n_phi);
			avg_mat.rhoZ_overA /= (double)(n_r*n_z*n_phi);
			avg_mat.rhoZ_overA_logI /= (double)(n_r*n_z*n_phi);
			
			MatTable[ir][iz] = avg_mat;
		}
		cout<<"\r Filling Material table ... "<<100.0*(double)ir/(double)Nr<<"%       ";cout.flush();
	}
	cout <<"Done"<<endl;
	
	table_initialized=true;
	pthread_mutex_unlock(&mutex);
}

//---------------------------------
// InitDRGeom
//---------------------------------
void DRootGeom::InitDRGeom(void)
{
	if(!gGeoManager)new TGeoManager();
	DRGeom = hddsroot();
}

//---------------------------------
// FindNode
//---------------------------------
TGeoNode* DRootGeom::FindNode(double *x)
{
	pthread_mutex_lock(&mutex);
	
	if(!DRGeom)InitDRGeom();

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
	if(!DRGeom)InitDRGeom();
  TGeoNode *cnode = DRGeom->FindNode(x[0],x[1],x[2]);
  TGeoVolume *cvol = cnode->GetVolume();
	pthread_mutex_unlock(&mutex);

  return cvol;  

}

//---------------------------------
// FindMat
//---------------------------------
jerror_t DRootGeom::FindMat(DVector3 pos, double &density, double &A, double &Z, double &RadLen) const
{
  return FindMatLL(pos, density, A, Z, RadLen);
}

//---------------------------------
// FindMat
//---------------------------------
jerror_t DRootGeom::FindMat(DVector3 pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const
{
	return FindMatTable(pos, rhoZ_overA, rhoZ_overA_logI, RadLen);
}

//---------------------------------
// FindMatTable
//---------------------------------
jerror_t DRootGeom::FindMatTable(DVector3 pos,double &density, double &A, double &Z, double &RadLen) const
{
	if(!table_initialized)((DRootGeom*)this)->InitTable();

	// For now, this just finds the bin in the material map the given position is in
	// (i.e. no interpolation )
	double r = pos.Perp();
	double z = pos.Z();
	int ir = (int)floor((r-r0)/dr);
	int iz = (int)floor((z-z0)/dz);
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
// FindMatTable
//---------------------------------
jerror_t DRootGeom::FindMatTable(DVector3 pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const
{
	if(!table_initialized)((DRootGeom*)this)->InitTable();

	// For now, this just finds the bin in the material map the given position is in
	// (i.e. no interpolation )
	double r = pos.Perp();
	double z = pos.Z();
	int ir = (int)floor((r-r0)/dr);
	int iz = (int)floor((z-z0)/dz);
	if(ir<0 || ir>=Nr || iz<0 || iz>=Nz){
		rhoZ_overA = rhoZ_overA_logI = RadLen = 0.0;
		return RESOURCE_UNAVAILABLE;
	}
	
	const VolMat &mat = MatTable[ir][iz];
	rhoZ_overA = mat.rhoZ_overA;
	rhoZ_overA_logI = mat.rhoZ_overA_logI;
	RadLen = mat.RadLen;

	return NOERROR;
}

// Find material properties by material name
jerror_t DRootGeom::FindMat(const char* matname,double &rhoZ_overA,
			    double &rhoZ_overA_logI, double &RadLen) const{  
  pthread_mutex_lock(const_cast<pthread_mutex_t*>(&mutex));
  
  if(!DRGeom)((DRootGeom*)this)->InitDRGeom(); // cast away constness to ensure DRGeom is set
  
  TGeoMaterial *mat=DRGeom->GetMaterial(matname);
  if (mat==NULL){
    _DBG_<<"Missing material " << matname <<endl;
    return RESOURCE_UNAVAILABLE;
  }
  double A=mat->GetA();
  double Z=mat->GetZ();
  double density=mat->GetDensity();
  rhoZ_overA=density*Z/A;
  RadLen=mat->GetRadLen();

  // Get mean excitation energy.  For a mixture this is calculated according 
  // to Leo 2nd ed p 29 eq 2.42
  double LnI=0.;
  if (mat->IsMixture()) {
    const TGeoMixture * mixt = dynamic_cast <const TGeoMixture*> (mat);
    for (int i=0;i<mixt->GetNelements();i++){
      double w_i=mixt->GetWmixt()[i];
      double Z_i=mixt->GetZmixt()[i];
      // Mean excitation energy for the element; see Leo 2nd ed., p25
      double I_i=7.0+12.0*Z_i; //eV
      if (Z_i>=13) I_i=9.76*Z_i+58.8*pow(Z_i,-0.19);
      I_i*=1e-9; // convert to GeV
      LnI+=w_i*Z_i*log(I_i)/Z;
    }
  }
  else{
    // mean excitation energy for the element
    double I=7.0+12.0*Z; //eV
    if (Z>=13) I=9.76*Z+58.8*pow(Z,-0.19);
    I*=1e-9; // Convert to GeV
    LnI=log(I);
  }
  rhoZ_overA_logI=rhoZ_overA*LnI;

  pthread_mutex_unlock(const_cast<pthread_mutex_t*>(&mutex));
  
  return NOERROR;
}


//---------------------------------
// FindMatLL
//---------------------------------
jerror_t DRootGeom::FindMatLL(DVector3 pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const
{
	/// This is a wrapper for the other FindMatLL method that returns density, A, and Z.
	/// It is here to make easy comparisons between the LL and table methods.

  double density, A, Z,LnI;
  FindMatLL(pos, density, A, Z, RadLen,LnI);
  rhoZ_overA = density*Z/A;
  rhoZ_overA_logI = rhoZ_overA*LnI;
  
  return NOERROR;
}

//---------------------------------
// FindMatLL
//---------------------------------
jerror_t DRootGeom::FindMatLL(DVector3 pos,double &density, double &A, double &Z,
			double &RadLen) const{
  density=RadLen=A=Z=0.;

	pthread_mutex_lock(const_cast<pthread_mutex_t*>(&mutex));
	
	if(!DRGeom)((DRootGeom*)this)->InitDRGeom(); // cast away constness to ensure DRGeom is set

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
// FindMatLL
//---------------------------------
jerror_t DRootGeom::FindMatLL(DVector3 pos,double &density, double &A, double &Z,
			      double &RadLen, double &LnI
			      ) const{
  density=RadLen=A=Z=LnI=0.;

	pthread_mutex_lock(const_cast<pthread_mutex_t*>(&mutex));
	
	if(!DRGeom)((DRootGeom*)this)->InitDRGeom(); // cast away constness to ensure DRGeom is set

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

  // Get mean excitation energy.  For a mixture this is calculated according 
  // to Leo 2nd ed p 29 eq 2.42
  if (cmat->IsMixture()) {
    const TGeoMixture * mixt = dynamic_cast <const TGeoMixture*> (cmat);
    LnI=0.;
    for (int i=0;i<mixt->GetNelements();i++){
      double w_i=mixt->GetWmixt()[i];
      double Z_i=mixt->GetZmixt()[i];
      // Mean excitation energy for the element; see Leo 2nd ed., p25
      double I_i=7.0+12.0*Z_i; //eV
      if (Z_i>=13) I_i=9.76*Z_i+58.8*pow(Z_i,-0.19);
      I_i*=1e-9; // convert to GeV
      LnI+=w_i*Z_i*log(I_i)/Z;
    }
  }
  else{
    // mean excitation energy for the element
    double I=7.0+12.0*Z; //eV
    if (Z>=13) I=9.76*Z+58.8*pow(Z,-0.19);
    I*=1e-9; // Convert to GeV
    LnI=log(I);
  }
 
  return NOERROR;
}



//---------------------------------
// FindMat
//---------------------------------
struct VolMat DRootGeom::FindMat(double *x)
{

	pthread_mutex_lock(&mutex);
	if(!DRGeom)InitDRGeom();
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


