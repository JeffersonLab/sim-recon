// $Id$
//
//    File: DMagneticFieldMapCalibDB.cc
// Created: Thu Jul 19 13:58:21 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp61.jlab.org 8.10.1 i386)
//

#include <cmath>
using namespace std;

#include "DMagneticFieldMapCalibDB.h"

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>


//---------------------------------
// DMagneticFieldMapCalibDB    (Constructor)
//---------------------------------
DMagneticFieldMapCalibDB::DMagneticFieldMapCalibDB(JApplication *japp)
{
	int runnumber = 1;
	jcalib = japp->GetJCalibration(runnumber);
	int Npoints = ReadMap(runnumber); 
	if(Npoints==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		japp->Quit();
	}
}

//---------------------------------
// DMagneticFieldMapCalibDB    (Constructor)
//---------------------------------
DMagneticFieldMapCalibDB::DMagneticFieldMapCalibDB(JCalibration *jcalib)
{
	this->jcalib = jcalib;
	if(ReadMap()==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		exit(1);
	} 
}

//---------------------------------
// ~DMagneticFieldMapCalibDB    (Destructor)
//---------------------------------
DMagneticFieldMapCalibDB::~DMagneticFieldMapCalibDB()
{

}

//---------------------------------
// ReadMap
//---------------------------------
int DMagneticFieldMapCalibDB::ReadMap(int runnumber, string context)
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
	
	string namepath = "Magnets/Solenoid/solenoid_1500";
	cout<<"Reading Magnetic field map from "<<namepath<<" ..."<<endl;
	vector< vector<float> > Bmap;
	jcalib->Get(namepath, Bmap);
	cout<<Bmap.size()<<" entries found (";
	
	// The map should be on a grid with equal spacing in x,y, and z.
	// Here we want to determine the number of points in each of these
	// dimensions and the range. Note that all of or maps right now are
	// 2-dimensional with extent only in x and z.
	// The easiest way to do this is to use a map<float, int> to make a
	// histogram of the entries by using the key to hold the extent
	// so that the number of entries will be equal to the number of
	// different values.
	map<float, int> xvals;
	map<float, int> yvals;
	map<float, int> zvals;
	xmin = ymin = zmin = 1.0E6;
	xmax = ymax = zmax = -1.0E6;
	for(unsigned int i=0; i<Bmap.size(); i++){
		vector<float> &a = Bmap[i];
		float &x = a[0];
		float &y = a[1];
		float &z = a[2];
		
		// Convert from inches to cm and shift in z by 26 inches to align
		// the origin with the GEANT defined one.
		x *= 2.54;
		y *= 2.54;
		z = (z-26.0)*2.54;
		
		xvals[x] = 1;
		yvals[y] = 1;
		zvals[z] = 1;
		if(x<xmin)xmin=x;
		if(y<ymin)ymin=y;
		if(z<zmin)zmin=z;
		if(x>xmax)xmax=x;
		if(y>ymax)ymax=y;
		if(z>zmax)zmax=z;
	}
	Nx = xvals.size();
	Ny = yvals.size();
	Nz = zvals.size();
	cout<<" Nx="<<Nx;
	cout<<" Ny="<<Ny;
	cout<<" Nz="<<Nz;
	cout<<" )"<<endl;
	
	// Create 3D vector so we can index the values by [x][y][z]
	vector<DBfieldPoint_t> zvec(Nz);
	vector< vector<DBfieldPoint_t> > yvec;
	for(int i=0; i<Ny; i++)yvec.push_back(zvec);
	for(int i=0; i<Nx; i++)Btable.push_back(yvec);
	
	// Distance between map points for r and z
	dx = (xmax-xmin)/(double)(Nx-1);
	dy = (ymax-ymin)/(double)(Ny-1);
	dz = (zmax-zmin)/(double)(Nz-1);

	// Copy values into Btable
	for(unsigned int i=0; i<Bmap.size(); i++){
		vector<float> &a = Bmap[i];
		int xindex = (int)floor((a[0]-xmin+dx/2.0)/dx); // the +dx/2.0 guarantees against round-off errors
		int yindex = (int)(Ny<2 ? 0:floor((a[1]-ymin+dy/2.0)/dy));
		int zindex = (int)floor((a[2]-zmin+dz/2.0)/dz);
		DBfieldPoint_t *b = &Btable[xindex][yindex][zindex];
		b->x = a[0];
		b->y = a[1];
		b->z = a[2];
		b->Bx = a[3];
		b->By = a[4];
		b->Bz = a[5];
	}
	
	// Calculate gradient at every point in the map.
	// The derivatives are calculated wrt to the *index* of the
	// field map, not the physical unit. This is because
	// the fractions used in interpolation are fractions
	// between indices.
	for(int index_x=0; index_x<Nx; index_x++){
		int index_x0 = index_x - (index_x>0 ? 1:0);
		int index_x1 = index_x + (index_x<Nx-1 ? 1:0);
		for(int index_z=1; index_z<Nz-1; index_z++){
			int index_z0 = index_z - (index_z>0 ? 1:0);
			int index_z1 = index_z + (index_z<Nz-1 ? 1:0);
			
			// Right now, all of our maps are 2-D with no y extent
			// so we comment out the following for-loop and hardwire
			// the index_y values to zero
			//for(int index_y=1; index_y<Ny-1; index_y++){
			//	int index_y0 = index_y-1;
			//	int index_y1 = index_y+1;
			if(Ny>1){
				_DBG_<<"Field map appears to be 3 dimensional. Code is currently"<<endl;
				_DBG_<<"unable to handle this. Exiting..."<<endl;
				exit(-1);
			}else{
				int index_y, index_y0, index_y1;
				index_y=index_y0=index_y1=0;

				DBfieldPoint_t *Bx0 = &Btable[index_x0][index_y][index_z];
				DBfieldPoint_t *Bx1 = &Btable[index_x1][index_y][index_z];
				DBfieldPoint_t *By0 = &Btable[index_x][index_y0][index_z];
				DBfieldPoint_t *By1 = &Btable[index_x][index_y1][index_z];
				DBfieldPoint_t *Bz0 = &Btable[index_x][index_y][index_z0];
				DBfieldPoint_t *Bz1 = &Btable[index_x][index_y][index_z1];
				
				DBfieldPoint_t *g = &Btable[index_x][index_y][index_z];
				g->dBxdx = (Bx1->Bx - Bx0->Bx)/(double)(index_x1-index_x0);
				g->dBxdy = (By1->Bx - By0->Bx)/(double)(index_y1-index_y0);
				g->dBxdz = (Bz1->Bx - Bz0->Bx)/(double)(index_z1-index_z0);

				g->dBydx = (Bx1->By - Bx0->By)/(double)(index_x1-index_x0);
				g->dBydy = (By1->By - By0->By)/(double)(index_y1-index_y0);
				g->dBydz = (Bz1->By - Bz0->By)/(double)(index_z1-index_z0);

				g->dBzdx = (Bx1->Bz - Bx0->Bz)/(double)(index_x1-index_x0);
				g->dBzdy = (By1->Bz - By0->Bz)/(double)(index_y1-index_y0);
				g->dBzdz = (Bz1->Bz - Bz0->Bz)/(double)(index_z1-index_z0);
			}
		}
	}

	return Bmap.size();
}

//---------------------------------
// GetField
//---------------------------------
void DMagneticFieldMapCalibDB::GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method) const
{
	/// This calculates the magnetic field at an arbitrary point
	/// in space using the field map read from the calibaration
	/// database. It interpolates between grid points using the
	/// gradient values calculated in ReadMap (called from the
	/// constructor).

	Bx = By = Bz = 0.0;

	if(Ny>1){
		_DBG_<<"Field map appears to be 3 dimensional. Code is currently"<<endl;
		_DBG_<<"unable to handle this. Treating as phi symmetric using y=0."<<endl;
	}

	// Get closest indices for this point
	double r = sqrt(x*x + y*y);
	int index_x = (int)floor((r-xmin)/dx + 0.5);
	if(index_x<0 || index_x>=Nx)return;
	int index_z = (int)floor((z-zmin)/dz + 0.5);	
	if(index_z<0 || index_z>=Nz)return;
	
	int index_y = 0;

	const DBfieldPoint_t *B = &Btable[index_x][index_y][index_z];

	// Fractional distance between map points.
	double ur = (r - B->x)/dx;
	double uz = (z - B->z)/dz;

	// Use gradient to project grid point to requested position
	double Br;
	Br = B->Bx + B->dBxdx*ur + B->dBxdz*uz;
	Bz = B->Bz + B->dBzdx*ur + B->dBzdz*uz;

	// Convert r back to x,y components
	double cos_theta = x/r;
	double sin_theta = y/r;
	if(r==0.0){
		cos_theta=1.0;
		sin_theta=0.0;
	}

	// Rotate back into phi direction
	Bx = Br*cos_theta;
	By = Br*sin_theta;
}


