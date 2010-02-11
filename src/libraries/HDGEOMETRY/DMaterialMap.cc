

#include <iostream>
using namespace std;

#include <DANA/DApplication.h>

#include "DMaterialMap.h"


//-----------------
// DMaterialMap  (Constructor)
//-----------------
DMaterialMap::DMaterialMap(string namepath, JCalibration *jcalib)
{
	/// Read the specified material map in from the calibration database.
	/// This will read in the map and figure out the number of grid
	/// points in each direction (r, and z) and the range in each.

	this->namepath = namepath;

	// Read in map from calibration database. This should really be
	// built into a run-dependent geometry framework, but for now
	// we do it this way.
	this->jcalib = jcalib;
	if(!jcalib)return;
	
	cout<<"Reading "<<namepath<<" ... "; cout.flush();
	vector< vector<float> > Mmap;
	jcalib->Get(namepath, Mmap);
	cout<<(int)Mmap.size()<<" entries (";
	if(Mmap.size()<1){
		cout<<")"<<endl;
		return;
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
	cout<<"Nr="<<Nr;
	cout<<" Nz="<<Nz;
	cout<<")"<<endl;
	
	// The values in the table are stored with r,z at the center of the node.
	// This means the actual map extends half a bin further out than the current
	// (local varfiable) rmin,rmax and zmin,zmax values. Set the class data
	// members to be the actual map limits
	this->rmin = rmin-dr/2.0;
	this->rmax = rmax+dr/2.0;
	this->zmin = zmin-dz/2.0;
	this->zmax = zmax+dz/2.0;
	
	// Set sizes of nested vectors to hold node data
	nodes.resize(Nr);
	for(int ir=0; ir<Nr; ir++){
		nodes[ir].resize(Nz);
	}
	
	// Fill table
	for(unsigned int i=0; i<Mmap.size(); i++){
		vector<float> &a = Mmap[i];
		float &r = a[0];
		float &z = a[1];
		int ir = (int)floor((r-this->rmin)/dr);
		int iz = (int)floor((z-this->zmin)/dz);
		if(ir<0 || ir>=Nr){_DBG_<<"ir out of range: ir="<<ir<<"  Nr="<<Nr<<endl; continue;}
		if(iz<0 || iz>=Nz){_DBG_<<"iz out of range: iz="<<iz<<"  Nz="<<Nz<<endl; continue;}
		MaterialNode &node = nodes[ir][iz];
		node.A = a[2];
		node.Z = a[3];
		node.Density = a[4];
		node.RadLen = a[5];
		node.rhoZ_overA = a[6];
		node.rhoZ_overA_logI = a[7];
	}
}


//-----------------
// FindNode
//-----------------
const DMaterialMap::MaterialNode* DMaterialMap::FindNode(DVector3 &pos) const
{
	// For now, this just finds the bin in the material map the given position is in
	// (i.e. no interpolation )
	double r = pos.Perp();
	double z = pos.Z();
	int ir = (int)floor((r-rmin)/dr);
	int iz = (int)floor((z-zmin)/dz);
	if(ir<0 || ir>=Nr || iz<0 || iz>=Nz)return NULL;
	
	return &nodes[ir][iz];
}

//-----------------
// FindMat
//-----------------
jerror_t DMaterialMap::FindMat(DVector3 &pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const
{
	const MaterialNode *node = FindNode(pos);
	if(!node)return RESOURCE_UNAVAILABLE;
	
	     rhoZ_overA = node->rhoZ_overA;
	rhoZ_overA_logI = node->rhoZ_overA_logI;
	         RadLen = node->RadLen;

	return NOERROR;
}

//-----------------
// FindMat
//-----------------
jerror_t DMaterialMap::FindMat(DVector3 &pos, double &density, double &A, double &Z, double &RadLen) const
{
	const MaterialNode *node = FindNode(pos);
	if(!node)return RESOURCE_UNAVAILABLE;
	
	density = node->Density;
			A = node->A;
			Z = node->Z;
	 RadLen = node->RadLen;

	return NOERROR;
}

//-----------------
// IsInMap
//-----------------
bool DMaterialMap::IsInMap(DVector3 &pos) const
{
	double r = pos.Perp();
	double z = pos.Z();
	return (r>=rmin) && (r<=rmax) && (z>=zmin) && (z<=zmax);
}
