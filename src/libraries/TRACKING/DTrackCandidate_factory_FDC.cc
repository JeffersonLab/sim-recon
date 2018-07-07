// $Id$
//
//    File: DTrackCandidate_factory.cc
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include <TROOT.h>

#include <math.h>

#include "DTrackCandidate_factory_FDC.h"
#include "DANA/DApplication.h"
#include "DQuickFit.h"
#include "JANA/JGeometry.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "DVector2.h"
#include "FDC/DFDCGeometry.h"
#include "DHoughFind.h"


bool FDCSortByZincreasing(DTrackCandidate_factory_FDC::DFDCTrkHit* const &hit1, DTrackCandidate_factory_FDC::DFDCTrkHit* const &hit2) {
	return hit1->hit->pos.Z() < hit2->hit->pos.Z();
}


//------------------
// DTrackCandidate_factory_FDC
//------------------
DTrackCandidate_factory_FDC::DTrackCandidate_factory_FDC()
{
#if 0
	// Set defaults
	MAX_SEED_DIST = 5.0;
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_DIST",			MAX_SEED_DIST);
#endif
}

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_FDC::init(void)
{
	TARGET_Z_MIN = 65.0 - 15.0;
	TARGET_Z_MAX = 65.0 + 15.0;
	
	MAX_HIT_DIST = 5.0;
	MAX_HIT_DIST2 = MAX_HIT_DIST*MAX_HIT_DIST;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_FDC::brun(JEventLoop *loop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_FDC::fini(void)
{	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_FDC::evnt(JEventLoop *loop, uint64_t eventnumber)
{

	// Get the hits into the trkhits vector
	GetTrkHits(loop);
	
	// Find seeds from X/Y projections
	vector<DFDCSeed> seeds;
	FindSeeds(seeds);
	
	// Loop over seeds and fit in phi-z plane to find z and theta
	for(unsigned int i=0; i<seeds.size(); i++){
		DFDCSeed &seed = seeds[i];
		if(debug_level>3)_DBG_<<"----- Seed "<<i<<" ------"<<endl;
		if(!seed.valid)continue;
		
		// Fit seed hits to get theta and vertex z position
		//FindThetaZ(seed);
		//if(!seed.valid)continue;

		// copy fit results to local variables (makes it easier to debug)
		double p_trans = seed.p_trans;
		double phi = seed.phi;
		double q = seed.q;
		double theta = seed.theta;
		double z_vertex = seed.z_vertex;

		if(debug_level>3)_DBG_<<"p_trans="<<p_trans<<" phi="<<phi<<" theta="<<theta<<" p="<<p_trans/sin(theta)<<endl;

		//Make a track candidate from results
		DTrackCandidate *can = new DTrackCandidate;
		DVector3 pos, mom;
		pos.SetXYZ(0.0, 0.0, z_vertex);
		mom.SetMagThetaPhi(p_trans/sin(theta), theta, phi);
		can->setPosition(pos);
		can->setMomentum(mom);
		can->setPID((q > 0.0) ? PiPlus : PiMinus);

		_data.push_back(can);
	}

	return NOERROR;
}

//------------------
// GetTrkHits
//------------------
void DTrackCandidate_factory_FDC::GetTrkHits(JEventLoop *loop)
{
	// Clear out old hits
	for(unsigned int i=0; i<fdctrkhits.size(); i++){
		delete fdctrkhits[i];
	}
	fdctrkhits.clear();

	// Get hits
	vector<const DFDCIntersection*> fdcintersects;
	loop->Get(fdcintersects);
	
	// Create a DFDCTrkHit object for each DFDCIntersection object
	for(unsigned int i=0; i<fdcintersects.size(); i++){
		
		DFDCTrkHit *trkhit = new DFDCTrkHit;
		trkhit->hit = fdcintersects[i];
		trkhit->flags = NOISE;
		
		fdctrkhits.push_back(trkhit);
	}
	
	// Filter out noise hits. All hits are initially flagged as "noise".
	// Hits with a neighbor within MAX_HIT_DIST have their noise flags cleared.
	for(unsigned int i=0; i<fdctrkhits.size(); i++){ // cut above should ensure cdctrkhits.size() is at least 1
		DFDCTrkHit *trkhit1 = fdctrkhits[i];
		for(unsigned int j=i+1; j<fdctrkhits.size(); j++){
			DFDCTrkHit *trkhit2 = fdctrkhits[j];
			if(!(trkhit2->flags & NOISE))continue; // this hit already has noise flag cleared
			double d2 = trkhit1->Dist2(trkhit2);
			double deltaz = fabs(trkhit1->hit->pos.Z() - trkhit2->hit->pos.Z());
			if(debug_level>4)_DBG_<<" -- Dist hits "<<i<<" and "<<j<<" : deltaR="<<sqrt(d2)<<"  deltaZ="<<deltaz<<endl;
			if(d2<MAX_HIT_DIST2 && deltaz<12.0){
				trkhit1->flags &= ~NOISE;
				trkhit2->flags &= ~NOISE;
			}// if
		}// j
	}// i
}



//------------------
// FindSeeds
//------------------
void DTrackCandidate_factory_FDC::FindSeeds(vector<DFDCSeed> &seeds)
{
	// Create a DHoughFind object to find a seed via the Hough Transform method
	// We create it once here to avoid the overhead of reallocating memory
	// each time through the loop below.
	//DHoughFind hough(-400.0, +400.0, -400.0, +400.0, 100, 100);

	// Loop until we find all seeds
	unsigned int last_available_hits = 0;
	for(int i=0; i<12; i++){ // limit the number of iterations in this loop
		// Check if we should exit the loop due to lack of available
		// hits. This is the primary point of exit for this loop.
		unsigned int Navailable_hits = NumAvailableHits();
		if(debug_level>3)_DBG_<<"Navailable_hits="<<Navailable_hits<<"  last_available_hits="<<last_available_hits<<endl;
		if(Navailable_hits<3)break;
		if(Navailable_hits==last_available_hits)break;
		last_available_hits = Navailable_hits;
	
		// Set the limits for the rough Hough histogram
		hough.SetLimits(-400.0, +400.0, -400.0, +400.0, 100, 100);
		hough.ClearPoints();

		// Add all points not already marked as unavailable
		for(unsigned int i=0; i<fdctrkhits.size(); i++){
			DFDCTrkHit *fdctrkhit = fdctrkhits[i];
			if(fdctrkhit->flags&USED)continue;
			if(fdctrkhit->flags&NOISE)continue;
			if(fdctrkhit->flags&OUT_OF_TIME)continue;
			const DFDCIntersection *hit = fdctrkhit->hit;
		
			hough.AddPoint(hit->pos.X(), hit->pos.Y());
		}
	
		DVector2 Ro = hough.Find();
		if(debug_level>3)_DBG_<<"Rough: GetMaxBinContent="<<hough.GetMaxBinContent()<<"  x0="<<Ro.X()<<"  y0="<<Ro.Y()<<endl;
		if(debug_level>6)hough.PrintHist();
		if(hough.GetMaxBinContent()<10.0)break;
		
		// Zoom in on resonanace a little
		double width = 60.0;
		hough.SetLimits(Ro.X()-width, Ro.X()+width, Ro.Y()-width, Ro.Y()+width, 100, 100);
		Ro = hough.Find();
		if(debug_level>3)_DBG_<<"Medium: GetMaxBinContent="<<hough.GetMaxBinContent()<<"  x0="<<Ro.X()<<"  y0="<<Ro.Y()<<endl;
		if(debug_level>5)hough.PrintHist();

		// Zoom in on resonanace once more
		width = 8.0;
		hough.SetLimits(Ro.X()-width, Ro.X()+width, Ro.Y()-width, Ro.Y()+width, 100, 100);
		Ro = hough.Find();
		if(debug_level>3)_DBG_<<"Fine: GetMaxBinContent="<<hough.GetMaxBinContent()<<"  x0="<<Ro.X()<<"  y0="<<Ro.Y()<<endl;
		if(debug_level>5)hough.PrintHist();
		
		DFDCSeed seed;
		
		seed.valid = true;
		seed.q = +1.0;
		seed.z_vertex = 65.0;
		seed.theta = M_PI/2.0;

		seed.r0 = Ro.Mod();
		seed.x0 = Ro.X();
		seed.y0 = Ro.Y();
		seed.p_trans = 0.003*2.1*seed.r0;
		
		// Set the VALID_HIT flag on hits consistent with this circle
		// and calculate the phi angle relative to the center of the
		// circle for the valid hits.
		FillSeedHits(seed);
		if(seed.hits.size()<3)continue; // require at least 5 hits

		// Set phi here since FillSeedHits may flip the sign of seed.q
		seed.phi = Ro.Phi() - seed.q*M_PI/2.0;

		// Find the theta and z-vertex of the seed and flag any used hits as USED
		FindThetaZ(seed); 

		// Add to list of seeds (i.e. tracks)
		seeds.push_back(seed);
	}
}

//------------------
// FillSeedHits
//------------------
void DTrackCandidate_factory_FDC::FillSeedHits(DFDCSeed &seed)
{
	// Loop over available hits and add any consistent with the circle
	// parameters to the list of hits for this seed.
	DVector2 Ro(seed.x0, seed.y0);
	seed.hits.clear();
	for(unsigned int i=0; i<fdctrkhits.size(); i++){
		DFDCTrkHit *fdctrkhit = fdctrkhits[i];
		fdctrkhit->flags &= ~ON_CIRCLE; // clear ON_CIRCLE flag, it gets set below if applicable
		if(fdctrkhit->flags&USED)continue;
		if(fdctrkhit->flags&NOISE)continue;
		if(fdctrkhit->flags&OUT_OF_TIME)continue;
		const DFDCIntersection *hit = fdctrkhit->hit;

		// Calculate distance between Hough transformed line (i.e.
		// the line on which a circle that passes through both the
		// origin and the point at hit->pos) and the circle center.
		DVector2 h(hit->pos.X()/2.0, hit->pos.Y()/2.0);
		DVector2 g(h.Y(), -h.X()); 
		g /= g.Mod();
		DVector2 Ro_minus_h(seed.x0-h.X(), seed.y0-h.Y());
		//Ro_minus_h /= Ro_minus_h.Mod();
		double dist = fabs(g.X()*Ro_minus_h.Y() - g.Y()*Ro_minus_h.X());

		// If this is not close enough to the found circle's center,
		// reject it for this seed.
		if(debug_level>3)_DBG_<<"dist="<<dist<<endl;
		if(dist > 2.0){
			fdctrkhit->flags &= ~VALID_HIT;
			continue;
		}
		
		fdctrkhit->flags |= ON_CIRCLE | VALID_HIT;
		seed.hits.push_back(fdctrkhit);
	}

	// Sort hits by increasing z
	stable_sort(seed.hits.begin(), seed.hits.end(), FDCSortByZincreasing);
	
	// Below, we fill in the phi_hit field of the DFDCTrkHit objects for
	// this candidate. We do so incrementally to try and accomodate
	// phi values greater than 2pi. Initialize the "last" direction as
	// pointing back to the beamline.
	DVector2 last_dir = -1.0*Ro/Ro.Mod();
	double last_phi = 0.0;

	// Loop over hits, filling in phi_hit
	int Nq = 0;
	for(unsigned int i=0; i<seed.hits.size(); i++){
		DFDCTrkHit *fdctrkhit = seed.hits[i];
		const DVector3 &pos = fdctrkhit->hit->pos;
		
		// Calculate phi. We do this trivially for now
		DVector2 p(pos.X() , pos.Y());
		DVector2 p_minus_Ro = p - Ro;
		p_minus_Ro/=p_minus_Ro.Mod();
		double delta_phi = p_minus_Ro.Phi() - last_dir.Phi();
		while(delta_phi>+M_PI)delta_phi-=2.0*M_PI;
		while(delta_phi<-M_PI)delta_phi+=2.0*M_PI;
		fdctrkhit->phi_hit = last_phi + delta_phi;
		last_phi = fdctrkhit->phi_hit;
		last_dir = p_minus_Ro;
		if(debug_level>3)_DBG_<<"phi_hit="<<fdctrkhit->phi_hit<<" z="<<fdctrkhit->hit->pos.Z()<<"  phi_hit/z="<<fdctrkhit->phi_hit/(fdctrkhit->hit->pos.Z()-65.0)<<" delta_phi="<<delta_phi<<endl;
		if(fdctrkhit->hit->wire2->layer<=6){
			Nq += fdctrkhit->phi_hit>0.0 ? +1:-1;
		}
	}
	
	// If Nq has too few entries, then look to the second package
	if(abs(Nq)<2){
		for(unsigned int i=0; i<seed.hits.size(); i++){
			DFDCTrkHit *fdctrkhit = seed.hits[i];
			if(fdctrkhit->hit->wire2->layer<=12){
				Nq += fdctrkhit->phi_hit>0.0 ? +1:-1;
			}
		}
	}
	
	if(Nq<0)seed.q = -seed.q;
}

//------------------
// NumAvailableHits
//------------------
unsigned int DTrackCandidate_factory_FDC::NumAvailableHits(void)
{
	// Loop over all hits and count the ones that are still
	// available for making a new seed.
	unsigned int N=0;
	for(unsigned int i=0; i<fdctrkhits.size(); i++){
		DFDCTrkHit *fdctrkhit = fdctrkhits[i];
		if(fdctrkhit->flags&USED)continue;
		if(fdctrkhit->flags&NOISE)continue;
		if(fdctrkhit->flags&OUT_OF_TIME)continue;
		if(fdctrkhit->flags&CANT_BE_IN_SEED)continue;
		N++;
	}
	
	return N;
}

//------------------
// FindThetaZ
//------------------
void DTrackCandidate_factory_FDC::FindThetaZ(DFDCSeed &seed)
{
	FindTheta(seed, TARGET_Z_MIN, TARGET_Z_MAX);
	FindZ(seed, seed.theta_min, seed.theta_max);
	
	// If z_vertex is not inside the target limits, then flag this
	// seed as invalid.
	if(seed.z_vertex<TARGET_Z_MIN || seed.z_vertex>TARGET_Z_MAX){
		if(debug_level>3)_DBG_<<"Seed z-vertex outside of target range (z="<<seed.z_vertex<<" TARGET_Z_MIN="<<TARGET_Z_MIN<<" TARGET_Z_MAX="<<TARGET_Z_MAX<<endl;
		seed.valid=false;
	}
	
	// If the seed is valid, mark all hits that are consistent with the track
	// as USED
	if(seed.valid){
		for(unsigned int i=0; i<seed.hits.size(); i++){
			DFDCTrkHit *trkhit = seed.hits[i];
			if(trkhit->flags&IN_Z_RANGE)trkhit->flags |= USED;
		}
	}else{
		// If the seed is not valid, then we need to mark at least one of the
		// hits us unusable on the next seed. Otherwise, the same seed will be found
		// over and over and over... We just mark the first hit in the seed.
		for(unsigned int i=0; i<seed.hits.size(); i++){
			DFDCTrkHit *trkhit = seed.hits[i];
			if(trkhit->flags&ON_CIRCLE)trkhit->flags |= CANT_BE_IN_SEED;
		}
	}
	
	return;
}

//------------------
// FindTheta
//------------------
void DTrackCandidate_factory_FDC::FindTheta(DFDCSeed &seed, double target_z_min, double target_z_max)
{
	/// Find the theta value using the hits from <i>seed</i>. 
	/// The value of seed.r0 is used to calculate theta.
	///
	/// This uses a histogramming technique that looks at the overlaps of the
	/// angle ranges subtended by each hit between the given target limits.
	/// The overlaps usually lead to a range of values for theta. The limits
	/// of these are stored in the theta_min and theta_max fields of the seed.
	/// The centroid of the range is stored in the theta field.
	
	// We use a simple array to store our histogram here. We don't want to use
	// ROOT histograms because they are not thread safe.
	unsigned int Nbins = 1000;
	unsigned int hist[Nbins];
	for(unsigned int i=0; i<Nbins; i++)hist[i] = 0; // clear histogram
	double bin_width = 2.0*M_PI/(double)Nbins;
	double hist_low_limit = -M_PI; // lower edge of histogram limits
	
	// Loop over valid hits, filling the histogram
	double &r0 = seed.r0;
	for(unsigned int i=0; i<seed.hits.size(); i++){
		DFDCTrkHit *trkhit = seed.hits[i];
		if(!trkhit->flags&VALID_HIT)continue;
		if(!trkhit->flags&ON_CIRCLE)continue;
		
		// Calculate upper and lower limits in theta
		double alpha = r0*trkhit->phi_hit;
		if(seed.q<0.0)alpha = -alpha;
		double z_hit = trkhit->hit->pos.Z();
		double tmin = atan2(alpha, z_hit - target_z_min);
		double tmax = atan2(alpha, z_hit - target_z_max);
		if(tmin>tmax){
			double tmp = tmin;
			tmin=tmax;
			tmax=tmp;
		}
		if(debug_level>3)_DBG_<<" -- phi_hit="<<trkhit->phi_hit<<" z_hit="<<z_hit<<endl;
		if(debug_level>3)_DBG_<<" -- tmin="<<tmin<<"  tmax="<<tmax<<endl;
		
		// Copy theta limits into trkhit
		trkhit->theta_min = tmin;
		trkhit->theta_max = tmax;
		
		// Find index of bins corresponding to tmin and tmax
		unsigned int imin = (unsigned int)floor((tmin-hist_low_limit)/bin_width);
		unsigned int imax = (unsigned int)floor((tmax-hist_low_limit)/bin_width);
		
		// If entire range of this hit is outside of the histogram limit
		// then discard this hit.
		if(imin>=Nbins)continue;
		
		// Clip limits of bin range to our histogram limits
		if(imin>=Nbins)imin=Nbins-1;
		if(imax>=Nbins)imax=Nbins-1;
		
		// Increment all bins between imin and imax
		for(unsigned int j=imin; j<=imax; j++)hist[j]++;
	}
	
	// Look for the indexes of the plateau
	unsigned int istart=0;
	unsigned int iend=0;
	for(unsigned int i=1; i<Nbins; i++){
		if(hist[i]>hist[istart]){
			istart = i;
			if(debug_level>3)_DBG_<<" -- istart="<<istart<<" (theta="<<hist_low_limit + bin_width*(0.5+(double)istart)<<" , N="<<hist[i]<<")"<<endl;
		}
		if(hist[i] == hist[istart])iend = i;
	}
	
	// If there are no entries in the histogram, then flag this seed as invalid
	if(hist[istart]==0.0){
		if(debug_level>3)_DBG_<<" - No entries in theta hist. Seed aborted."<<endl;
		seed.valid=false;
	}
	
	// Calculate theta limits
	seed.theta_min = hist_low_limit + bin_width*(0.5+(double)istart);
	seed.theta_max = hist_low_limit + bin_width*(0.5+(double)iend);
	seed.theta = (seed.theta_max + seed.theta_min)/2.0;
	if(debug_level>3)_DBG_<<"istart="<<istart<<" iend="<<iend<<" theta_min="<<seed.theta_min<<" theta_max="<<seed.theta_max<<endl;

	// Mark all of the on circle hits that have a theta range consistent with seed.theta
	// as being IN_THETA_RANGE.
	for(unsigned int i=0; i<seed.hits.size(); i++){
		DFDCTrkHit *trkhit = seed.hits[i];
		trkhit->flags &= ~IN_THETA_RANGE;
		if(!trkhit->flags&VALID_HIT)continue;
		if(!trkhit->flags&ON_CIRCLE)continue;
		if(trkhit->theta_min > seed.theta)continue;
		if(trkhit->theta_max < seed.theta)continue;
		trkhit->flags |= IN_THETA_RANGE;
	}

	// If theta is negative, then we probably chose the wrong sign. Flip
	// it now.
	if(seed.theta<0.0){
		seed.theta_min *= -1.0;
		seed.theta_max *= -1.0;
		seed.theta *= -1.0;
		seed.q *= -1.0;
		seed.phi += M_PI;
		if(seed.phi > 2.0*M_PI)seed.phi -= 2.0*M_PI;
	}
}

//------------------
// FindZ
//------------------
void DTrackCandidate_factory_FDC::FindZ(DFDCSeed &seed, double theta_min, double theta_max)
{
	/// Find the z value of the vertex using the valid stereo hits from <i>seed</i>. The values
	/// for phi_hit and pos.Z() are assumed to be valid as is the status of the
	/// VALID_HIT bit in flags.
	///
	/// This uses a histogramming technique that looks at the overlaps of the
	/// z ranges subtended by each hit between the given theta limits.
	/// The overlaps usually lead to a range of values for z_vertex. The limits
	/// of these are stored in the z_min and z_max fields of the seed.
	/// The centroid of the range is stored in the z_vertex field.
	
	// We use a simple array to store our histogram here. We don't want to use
	// ROOT histograms because they are not thread safe.
	unsigned int Nbins = 300;
	unsigned int hist[Nbins];
	for(unsigned int i=0; i<Nbins; i++)hist[i] = 0; // clear histogram
	double bin_width = 0.5; // bins are 5mm
	double hist_low_limit = 0.0; // lower edge of histogram limits

	// We effectively extend the theta_min and theta_max angles here
	// a bit to include some error. The motivation is that if
	// theta_min == theta_max that leads to z_min == z_max so there
	// is little or no overlap of the z ranges of separate hits.
	// For now, we hardwire this to 1 degree
	double theta_err = 1.0/57.3;
	
	// Loop over valid hits, filling the histogram
	double r0 = seed.r0;
	double tan_alpha_max = tan(theta_max + theta_err)/r0;
	double tan_alpha_min = tan(theta_min - theta_err)/r0;
	if(tan_alpha_min<0.0)tan_alpha_min=0.0;
	for(unsigned int i=0; i<seed.hits.size(); i++){
		DFDCTrkHit *trkhit = seed.hits[i];
		if(!trkhit->flags&VALID_HIT)continue;
		if(!trkhit->flags&ON_CIRCLE)continue;
		if(!trkhit->flags&IN_THETA_RANGE)continue;
		
		// Calculate upper and lower limits in z
		double q_sign = seed.q>0.0 ? +1.0:-1.0;
		double z_hit = trkhit->hit->pos.Z();
		double zmin = z_hit - q_sign*trkhit->phi_hit/tan_alpha_min;
		double zmax = z_hit - q_sign*trkhit->phi_hit/tan_alpha_max;
		if(zmin>zmax){
			double tmp = zmin;
			zmin=zmax;
			zmax=tmp;
		}
		
		// Copy z limits into trkhit
		trkhit->zmin = zmin;
		trkhit->zmax = zmax;
		
		// Find index of bins corresponding to tmin and tmax
		unsigned int imin = zmin<hist_low_limit ? 0:(unsigned int)floor((zmin-hist_low_limit)/bin_width);
		unsigned int imax = zmax<hist_low_limit ? 0:(unsigned int)floor((zmax-hist_low_limit)/bin_width);
		
		if(debug_level>3)_DBG_<<" -- phi_hit="<<trkhit->phi_hit<<" z_hit="<<z_hit<<endl;
		if(debug_level>3)_DBG_<<" -- zmin="<<zmin<<"  zmax="<<zmax<<endl;
		if(debug_level>3)_DBG_<<" -- imin="<<imin<<"  imax="<<imax<<endl;

		// If entire range of this hit is outside of the histogram limit
		// then discard this hit.
		if(imax<=0 || imin>=Nbins)continue;
		
		// Clip limits of bin range to our histogram limits
		if(imin>=Nbins)imin=Nbins-1;
		if(imax>=Nbins)imax=Nbins-1;
		
		// Increment all bins between imin and imax
		for(unsigned int j=imin; j<=imax; j++)hist[j]++;
	}
	
	// Look for the indexes of the plateau
	unsigned int istart=(unsigned int)((TARGET_Z_MIN-hist_low_limit)/bin_width);
	unsigned int iend=0;
	for(unsigned int i=1; i<Nbins; i++){
		
		// Only look in Target area
		double z = hist_low_limit + bin_width*(0.5+(double)i);
		if(z<TARGET_Z_MIN || z>TARGET_Z_MAX)continue;
	
		if(hist[i]>hist[istart]){
			istart = i;
			if(debug_level>3)_DBG_<<" -- istart="<<istart<<" (z="<<hist_low_limit + bin_width*(0.5+(double)istart)<<" , N="<<hist[i]<<")"<<endl;
		}
		if(hist[i] == hist[istart])iend = i;
	}

	// If there are no entries in the histogram, then flag this seed as invalid
	if(hist[istart]==0.0){
		if(debug_level>3)_DBG_<<" - No entries in z-vertex hist. Seed aborted."<<endl;
		seed.valid=false;
	}
	
	// Calculate z limits
	seed.z_min = hist_low_limit + bin_width*(0.5+(double)istart);
	seed.z_max = hist_low_limit + bin_width*(0.5+(double)iend);
	seed.z_vertex = (seed.z_max + seed.z_min)/2.0;
	if(debug_level>3)_DBG_<<"istart="<<istart<<" iend="<<iend<<" z_min="<<seed.z_min<<" z_max="<<seed.z_max<<" hits[istart]="<<hist[istart]<<endl;

	// Mark all of the hits that have a theta range consistent with seed.theta
	// and a z_vertex consistent with seed.z_vertex as being IN_Z_RANGE.
	for(unsigned int i=0; i<seed.hits.size(); i++){
		DFDCTrkHit *trkhit = seed.hits[i];
		trkhit->flags &= ~IN_Z_RANGE;
		if(!trkhit->flags&VALID_HIT)continue;
		if(!trkhit->flags&ON_CIRCLE)continue;
		if(!trkhit->flags&IN_THETA_RANGE)continue;
		if(trkhit->zmin > seed.z_vertex)continue;
		if(trkhit->zmax < seed.z_vertex)continue;
		trkhit->flags |= IN_Z_RANGE;
	}
}

