

#include <iostream>
#include <set>
using namespace std;

#include <DANA/DApplication.h>
#include <DVector2.h>

#include "DMaterialMap.h"
#include "DMagneticFieldMap.h"


//-----------------
// DMaterialMap  (Constructor)
//-----------------
DMaterialMap::DMaterialMap(string namepath, JCalibration *jcalib)
{
	/// Read the specified material map in from the calibration database.
	/// This will read in the map and figure out the number of grid
	/// points in each direction (r, and z) and the range in each.

	MAX_BOUNDARY_SEARCH_STEPS = 30;
	ENABLE_BOUNDARY_CHECK = true;
	
	gPARMS->SetDefaultParameter("GEOM:MAX_BOUNDARY_SEARCH_STEPS", MAX_BOUNDARY_SEARCH_STEPS, "Maximum number of steps (cells) to iterate when searching for a material boundary in DMaterialMap::EstimatedDistanceToBoundary(...)");
	gPARMS->SetDefaultParameter("GEOM:ENABLE_BOUNDARY_CHECK", ENABLE_BOUNDARY_CHECK, "Enable boundary checking (superceeds any setting in DReferenceTrajectory). This is for debugging only.");

	this->namepath = namepath;

	// Read in map from calibration database. This should really be
	// built into a run-dependent geometry framework, but for now
	// we do it this way.
	this->jcalib = jcalib;
	if(!jcalib)return;
	
	string blank(' ',80);
	cout<<blank<<"\r"; cout.flush(); // clear line
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
	// (local variable) rmin,rmax and zmin,zmax values. Set the class data
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
		node.LogI=node.rhoZ_overA_logI/(node.rhoZ_overA==0.0 ? 1.0:node.rhoZ_overA);
		node.KrhoZ_overA=0.1535e-3*node.rhoZ_overA;
	}
	
	// Now find the r and z boundaries that will be used during swimming
	FindBoundaries();
}

//-----------------
// FindBoundaries
//-----------------
void DMaterialMap::FindBoundaries(void)
{
	/// Look for boundaries where the material density changes
	/// and record them for faster checking during swimming
	
	// Boundaries are recorded for r(or z) bins for which ANY adjacent
	// z(or r) bin shows a factor of 2 or more change in the 
	// radiation length of the material.
	//
	// This is done by making 2 passes over every bin in the map
	// (one for the r boundaries and one for the z) and recording
	// the index of the bin for which an increase was detected.
	// Doing it this way will detect a boundary multiple times
	// (one for each bin in the direction the boundary runs).
	// Thus, we use a set container to record the boundaries
	// temporarily since writing the same value to it multiple
	// times won't result in multiple entries.
	//
	// A final pass at the end records the r (or z) boundary values
	// for the border between the bins in the vectors containers
	// that are actually used by the EstimatedDistanceToBoundary
	// method below.
	
	// Loop over r bins
	set<int> iz_boundaries;
	for(int ir=0; ir<Nr; ir++){
		
		// initialize boundary
		double RadLen1 = nodes[ir][0].RadLen;
		
		// Loop over z bins
		for(int iz=1; iz<Nz; iz++){
			double RadLen2 = nodes[ir][iz].RadLen;
			if(RadLen1<0.5*RadLen2 || RadLen2<0.5*RadLen1){
				iz_boundaries.insert(iz);
			}
		}
	}

	// Loop over z bins
	set<int> ir_boundaries;
	for(int iz=0; iz<Nz; iz++){
		
		// initialize boundary
		double RadLen1 = nodes[0][iz].RadLen;
		
		// Loop over r bins
		for(int ir=1; ir<Nr; ir++){
			double RadLen2 = nodes[ir][iz].RadLen;
			if(RadLen1<0.5*RadLen2 || RadLen2<0.5*RadLen1){
				ir_boundaries.insert(ir);
			}
		}
	}

	// Fill r_boundaries and z_boundaries with actual values of
	// borders between the bins
	for(set<int>::iterator iter=iz_boundaries.begin(); iter!=iz_boundaries.end(); iter++){
		int iz = *iter;
		z_boundaries.push_back(zmin+dz*(double)iz);
	}
	for(set<int>::iterator iter=ir_boundaries.begin(); iter!=ir_boundaries.end(); iter++){
		int ir = *iter;
		r_boundaries.push_back(rmin+dr*(double)ir);
	}


	// If the number of boundaries in either direction is more
	// than half the number of bins then flag this material
	// as having an irregular density profile.
	irregular_density_profile = false;
	irregular_density_profile |= (int)r_boundaries.size() > (Nr/2);
	irregular_density_profile |= (int)z_boundaries.size() > (Nz/2);
	if(irregular_density_profile){
		r_boundaries.clear();
		z_boundaries.clear();
	}
}

//-----------------
// FindMat
//-----------------
jerror_t DMaterialMap::FindMat(DVector3 &pos,double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const
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
jerror_t DMaterialMap::FindMatALT1(DVector3 &pos,double &KrhoZ_overA,
			       double &rhoZ_overA, double &logI, 
			       double &RadLen) const
{
	const MaterialNode *node = FindNode(pos);
	if(!node)return RESOURCE_UNAVAILABLE;
	
	KrhoZ_overA=node->KrhoZ_overA;
	rhoZ_overA = node->rhoZ_overA;
	logI = node->LogI;
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

// The Kalman filter needs a slightly different variant of the material 
// properties
jerror_t DMaterialMap::FindMatKalman(DVector3 &pos,double &Z,
				     double &K_rho_Z_over_A,
				     double &rho_Z_over_A,
				     double &LogI) const
{
	const MaterialNode *node = FindNode(pos);
	if(!node)return RESOURCE_UNAVAILABLE;
	
	Z = node->Z;
	rho_Z_over_A = node->rhoZ_overA;
	LogI = node->LogI;
	K_rho_Z_over_A=node->KrhoZ_overA;

	return NOERROR;
}


//-----------------
// IsInMap
//-----------------
bool DMaterialMap::IsInMap(const DVector3 &pos) const
{
	double pos_x = pos.X();
	double pos_y = pos.Y();
	double r = sqrt(pos_x*pos_x + pos_y*pos_y);
	double z = pos.Z();
	return (r>=rmin) && (r<=rmax) && (z>=zmin) && (z<=zmax);
}

//-----------------
// EstimatedDistanceToBoundary
//-----------------
double DMaterialMap::EstimatedDistanceToBoundary(const DVector3 &pos, const DVector3 &mom)
{
	/// Give a very rough estimate of the distance a track at this position/momentum
	/// will travel before seeing either a significant change in the raditation
	/// length of material, or the edge of this map. The primary purpose of
	/// this is to help determine the appropriate step size when swimming 
	/// charged tracks. As we approach a boundary of materials, we want to take
	/// smaller steps to ensure the material is integrated properly.
	///
	/// This method can be called for points either inside or outside of this map.
	/// It can be that this map overlaps another so for points outside, we look
	/// for the point where the track would enter this map.
	///
	/// If a problem is encountered (e.g. point is outside of map and not pointing
	/// toward it) we return a value of 1.0E6

	// The method here is to use either dr/dz or dz/dr to project a straight line
	// in r/z space and find bins along that line that can be be checked for
	// increases in material density. This will clearly give a worse answer than
	// projecting a helix, but should be pretty quick and robust. Speed is critical
	// here since this is called for every step during swimming.

	double s_to_boundary = 1.0E6;
	if(!ENABLE_BOUNDARY_CHECK)return s_to_boundary; // low-level, catch-all opportunity to NOT do this.
	
	// This gets called often enough it is worth trying to bail as quickly
	// as possible. We check if the point is outside the boundary and headed
	// away from the map first.
	double z = pos.Z();
	double pz = mom.Z();
	if(pz>0.0)
		if(z>zmax)return s_to_boundary;
	else
		if(z<zmin)return s_to_boundary;

	double pos_x = pos.X();
	double pos_y = pos.Y();
	double mom_x = mom.X();
	double mom_y = mom.Y();
	double x_dot_p = pos_x*mom_x + pos_y*mom_y;
	double r = sqrt(pos_x*pos_x + pos_y*pos_y);
	if(x_dot_p>0.0)
		if(r>rmax)return s_to_boundary;
	else
		if(r<rmin)return s_to_boundary;
	

	double pr = sqrt(mom_x*mom_x + mom_y*mom_y) * (x_dot_p>0 ? +1.0:-1.0); // sign of pr depends on x_dot_p
	
	// Unit vector pointing in momentum direction in r-z space
	double mod = sqrt(pr*pr + pz*pz);
	if(mod<1.0E-6)return s_to_boundary; // for when momentum is purely in phi direction
	double p_hatR = pr/mod;
	double p_hatZ = pz/mod;
	DVector2 p_hat(p_hatR, p_hatZ);

	// Get shortest distance to boundary of entire map
	s_to_boundary = DistanceToBox(r, z, p_hatR, p_hatZ, rmin, rmax, zmin, zmax);

	// If point is outside of map, then return now. s_to_boundary has valid answer
	if(r<=rmin || r>=rmax || z<=zmin || z>=zmax)return s_to_boundary;
	
	// At this point we know that the point is inside of this map. We need to
	// search for boundaries within this map to find the closest one to the
	// track. If the map is marked as "irregular" then we do a bin by bin 
	// search for changes in density. Otherwise, we loop over the boundary
	// map already generated during construction.
	if(irregular_density_profile){
		return EstimatedDistanceToBoundarySearch(r, z, p_hatR, p_hatZ, s_to_boundary);
	}

	
	// Loop over boundaries stored in z_boundaries
	for(unsigned int i=0; i<z_boundaries.size(); i++){
		double dist = (z_boundaries[i]-z)/p_hatZ;
		if(!finite(dist))continue;
		if(dist<0.0)continue; // boundary is behind us
		if(dist>s_to_boundary)continue; // boundary intersection is already further than closest boundary
		
		// Now we must calculate whether this intersection point is
		// still inside the map.
		double my_r = r + dist*p_hatR;
		if(my_r>=rmin && my_r<=rmax)s_to_boundary = dist;
	}
	
	// Loop over boundaries stored in r_boundaries
	for(unsigned int i=0; i<r_boundaries.size(); i++){
		double dist = (r_boundaries[i]-r)/p_hatR;
		if(!finite(dist))continue;
		if(dist<0.0)continue; // boundary is behind us
		if(dist>s_to_boundary)continue; // boundary intersection is already further than closest boundary
		
		// Now we must calculate whether this intersection point is
		// still inside the map.
		double my_z = z + dist*p_hatZ;
		if(my_z>=zmin && my_z<=zmax)s_to_boundary = dist;
	}

	return s_to_boundary;
}

//-----------------
// EstimatedDistanceToBoundarySearch
//-----------------
double DMaterialMap::EstimatedDistanceToBoundarySearch(double r, double z, double p_hatR, double p_hatZ, double &s_to_boundary)
{
	/// Estimate the distance to the nearest boundary (marked by a
	/// a change in radiation length by a factor of 2 or more) by
	/// doing an exhaustive bin-by-bin search in the direction of
	/// the momentum vector in r-z space. This should only get called
	/// for maps with an irregular density profile (i.e. one whose
	/// boundaries do not appear to be in pure r or pure z directions.

	// If we got this far then the point is inside of this map. We need to check
	// if there is a jump in the radiation length of a cell that is in the path
	// of the trajectory so we can return that. Otherwise, we'll return the distance
	// passed in via s_to_boundary.

	// We want to form a vector pointing in the direction of momentum in r,z space
	// but with a magnitude such that we are taking a step across a single grid
	// point in either r or z, whichever is smallest.
	double scale_r=fabs(dr/p_hatR);
	double scale_z=fabs(dz/p_hatZ);
	double scale=1000.0;
	if(finite(scale_r) && scale_r<scale)scale = scale_r;
	if(finite(scale_z) && scale_z<scale)scale = scale_z;
	DVector2 delta_rz(scale*p_hatR, scale*p_hatZ);

	// Find radiation length of our starting cell
	int ir_start = (int)floor((r-rmin)/dr);
	int iz_start = (int)floor((z-zmin)/dz);
	double RadLen_start = nodes[ir_start][iz_start].RadLen;
	
	// Loop until we find a change of radiation length within this map or
	// until we hit the edge of our boundaries.
	DVector2 rzpos_start(r, z);
	DVector2 last_rzpos = rzpos_start;
	int last_ir = ir_start;
	int last_iz = iz_start;
	for(int Nsteps=0; Nsteps<MAX_BOUNDARY_SEARCH_STEPS; Nsteps++){ // limit us to looking only 10 grid points away for speed
		// Step to next cell
		DVector2 rzpos = last_rzpos + delta_rz;

		// Find indexes for this cell
		int ir = (int)floor((rzpos.X()-rmin)/dr);
		int iz = (int)floor((rzpos.Y()-zmin)/dz);

		// Check if we hit the boundary of the map and if so, simply return
		// the value found above for the distance to the map boundary.
		if(ir<0 || ir>=Nr || iz<0 || iz>=Nz)return s_to_boundary;

		// Check radiation length against start point's
		double RadLen = nodes[ir][iz].RadLen;
		if(RadLen <= 0.5*RadLen_start){
			double rmin_cell = (double)last_ir*dr + rmin;
			double rmax_cell = rmin_cell + dr;
			double zmin_cell = (double)last_iz*dz + zmin;
			double zmax_cell = zmin_cell + dz;
			double last_rzposX = last_rzpos.X();
			double last_rzposY = last_rzpos.Y();
			double s_to_cell = DistanceToBox(last_rzposX, last_rzposY, p_hatR, p_hatZ, rmin_cell, rmax_cell, zmin_cell, zmax_cell);
			if(s_to_cell==1.0E6)s_to_cell = 0.0;
			double my_s_to_boundary = (last_rzpos-rzpos_start).Mod() + s_to_cell;
			if(my_s_to_boundary < s_to_boundary)s_to_boundary = my_s_to_boundary;
			
			return s_to_boundary;
		}
		
		// Need to take another step. Update "last" variables
		last_iz = iz;
		last_ir = ir;
		last_rzpos = rzpos;
	}
	
	return s_to_boundary;
}

//-----------------
// DistanceToBox
//-----------------
double DMaterialMap::DistanceToBox(double &posx, double &posy, double &xdir, double &ydir, double xmin, double xmax, double ymin, double ymax)
{
	/// Given a point in 2-D space, a direction, and the limits of a box, find the
	/// closest distance to box edge in the given direction.
	
	// Calculate intersection distances to all 4 boundaries of map
	double dist_x1 = (xmin-posx)/xdir;
	double dist_x2 = (xmax-posx)/xdir;
	double dist_y1 = (ymin-posy)/ydir;
	double dist_y2 = (ymax-posy)/ydir;

	// Make list of all positive, finite distances that are
	// on border of box.
	double shortestDist = 1.0E6;
	if(finite(dist_x1) && dist_x1>=0.0 && dist_x1 < shortestDist){
		//double y = (pos + dist_x1*dir).Y();
		double y = posy + dist_x1*ydir;
		if(y>=ymin && y<=ymax) shortestDist = dist_x1;
	}
	if(finite(dist_x2) && dist_x2>=0.0 && dist_x2 < shortestDist){
		//double y = (pos + dist_x2*dir).Y();
		double y = posy + dist_x2*ydir;
		if(y>=ymin && y<=ymax) shortestDist = dist_x2;
	}
	if(finite(dist_y1) && dist_y1>=0.0 && dist_y1 < shortestDist){
		//double x = (pos + dist_y1*dir).X();
		double x = posx + dist_y1*xdir;
		if(x>=xmin && x<=xmax) shortestDist = dist_y1;
	}
	if(finite(dist_y2) && dist_y2>=0.0 && dist_y2 < shortestDist){
		//double x = (pos + dist_y2*dir).X();
		double x = posx + dist_y2*xdir;
		if(x>=xmin && x<=xmax) shortestDist = dist_y2;
	}
  	
	// Return shortest distance
	return shortestDist;
}

