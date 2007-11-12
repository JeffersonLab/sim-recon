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
#include "DTrack.h"
#include "DQuickFit.h"
#include "JANA/JGeometry.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "DVector2.h"
#include "FDC/DFDCGeometry.h"
#include "DHoughFind.h"

#if 0
bool FDCSortByRdecreasing(DTrackCandidate_factory_FDC::DFDCTrkHit* const &hit1, DTrackCandidate_factory_FDC::DFDCTrkHit* const &hit2) {
	double r1_2 = hit1->hit->x*hit1->hit->x + hit1->hit->y*hit1->hit->y;
	double r2_2 = hit2->hit->x*hit2->hit->x + hit2->hit->y*hit2->hit->y;

	return r1_2 > r2_2;
}
bool FDCSortByRincreasing(DTrackCandidate_factory_FDC::DFDCTrkHit* const &hit1, DTrackCandidate_factory_FDC::DFDCTrkHit* const &hit2) {
	return !FDCSortByRdecreasing(hit1, hit2);
}
bool FDCSortByZdecreasing(DTrackCandidate_factory_FDC::DFDCTrkHit* const &hit1, DTrackCandidate_factory_FDC::DFDCTrkHit* const &hit2) {
	return hit1->hit->wire->layer > hit2->hit->wire->layer;
}
#endif

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
	MAX_SEED_HITS = 10;
	MIN_SEED_HITS = 4;
	MAX_CIRCLE_DIST = 2.0;
	MAX_PHI_Z_DIST = 10.0;
	MIN_PHI_Z_HITS = 4;
	MAX_DEBUG_BUFFERS = 0;
	TARGET_Z_MIN = 50.0;
	TARGET_Z_MAX = 80.0;
	TRACKHIT_SOURCE = "MC";
	XY_NOISE_CUT = 2.0;
	MIN_HIT_Z = -100.0;
	MAX_HIT_Z = +360.0;
	EXCLUDE_STEREO = true;
	MIN_CANDIDATE_HITS = 6;
	DEBUG_HISTS = false;
	
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_DIST",			MAX_SEED_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_HITS",			MAX_SEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_SEED_HITS",			MIN_SEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_CIRCLE_DIST",		MAX_CIRCLE_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_PHI_Z_DIST",			MAX_PHI_Z_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_PHI_Z_HITS",			MIN_PHI_Z_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_DEBUG_BUFFERS",		MAX_DEBUG_BUFFERS);
	gPARMS->SetDefaultParameter("TRKFIND:TARGET_Z_MIN",			TARGET_Z_MIN);
	gPARMS->SetDefaultParameter("TRKFIND:TARGET_Z_MAX",			TARGET_Z_MAX);
	gPARMS->SetDefaultParameter("TRKFIND:TRACKHIT_SOURCE",		TRACKHIT_SOURCE);
	gPARMS->SetDefaultParameter("TRKFIND:XY_NOISE_CUT",			XY_NOISE_CUT);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_HIT_Z",				MIN_HIT_Z);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_HIT_Z",				MAX_HIT_Z);
	gPARMS->SetDefaultParameter("TRKFIND:EXCLUDE_STEREO",			EXCLUDE_STEREO);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_CANDIDATE_HITS",	MIN_CANDIDATE_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:DEBUG_HISTS",				DEBUG_HISTS);

	MAX_SEED_DIST2 = MAX_SEED_DIST*MAX_SEED_DIST;
	XY_NOISE_CUT2 = XY_NOISE_CUT*XY_NOISE_CUT;
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
jerror_t DTrackCandidate_factory_FDC::brun(JEventLoop *loop, int runnumber)
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
jerror_t DTrackCandidate_factory_FDC::evnt(JEventLoop *loop, int eventnumber)
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
		FindThetaZ(seed);
		if(!seed.valid)continue;

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
		can->setCharge(q);

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
		if(!(trkhit1->flags & NOISE))continue; // this hit already not marked for noise
		for(unsigned int j=0; j<fdctrkhits.size(); j++){
			if(j==i)continue;
			double d2 = trkhit1->Dist2(fdctrkhits[j]);
			if(d2<MAX_HIT_DIST2){
				trkhit1->flags &= ~NOISE;
				fdctrkhits[j]->flags &= ~NOISE;
				break;
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
	DHoughFind hough(-400.0, +400.0, -400.0, +400.0, 100, 100);

	// Loop until we find all seeds
	unsigned int last_available_hits = 0;
	do{
		// Check if we should exit the loop due to lack of available
		// hits. This is the primary point of exit for this loop.
		unsigned int Navailable_hits = NumAvailableHits();
		if(debug_level>3)_DBG_<<"Navailable_hits="<<Navailable_hits<<"  last_available_hits="<<last_available_hits<<endl;
		if(Navailable_hits<3)break;
		if(Navailable_hits==last_available_hits)break;
		last_available_hits = Navailable_hits;
	
		// Set the limits for the rough Hough histogram
		hough.SetLimits(-400.0, +400.0, -400.0, +400.0, 100, 100);

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
		if(hough.GetMaxBinContent()<10.0)continue;
		
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
		seed.phi = Ro.Phi() - seed.q*M_PI/2.0;
		seed.p_trans = 0.003*2.1*seed.r0;
		
		// Set the VALID_HIT flag on hits consistent with this circle
		// and calculate the phi angle relative to the center of the
		// circle for the vaid hits.
		FillSeedHits(seed);
		if(seed.hits.size()<3)continue; // require at least 5 hits
		
		seeds.push_back(seed);
	}while(true); // only loop once for now while debugging
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
		
		fdctrkhit->flags |= USED | VALID_HIT;
		seed.hits.push_back(fdctrkhit);
	}

	// Sort hits by increasing z
	sort(seed.hits.begin(), seed.hits.end(), FDCSortByZincreasing);
	
	// Below, we fill in the phi_hit field of the DFDCTrkHit objects for
	// this candidate. We do so incrementally to try and accomodate
	// phi values greater than 2pi. Initialize the "last" direction as
	// pointing back to the beamline.
	DVector2 last_dir = -1.0*Ro/Ro.Mod();
	double last_phi = 0.0;

	// Loop over hits, filling in phi_hit
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
	}
	
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
		
		// Find index of bins corresponding to tmin and tmax
		unsigned int imin = (unsigned int)floor((tmin-hist_low_limit)/bin_width);
		unsigned int imax = (unsigned int)floor((tmax-hist_low_limit)/bin_width);
		
		// If entire range of this hit is outside of the histogram limit
		// then discard this hit.
		if(imax<0 || imin>=Nbins)continue;
		
		// Clip limits of bin range to our histogram limits
		if(imin<0)imin=0;
		if(imin>=Nbins)imin=Nbins-1;
		if(imax<0)imax=0;
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
	if(hist[istart]==0.0)seed.valid=false;
	
	// Calculate theta limits
	seed.theta_min = hist_low_limit + bin_width*(0.5+(double)istart);
	seed.theta_max = hist_low_limit + bin_width*(0.5+(double)iend);
	seed.theta = (seed.theta_max + seed.theta_min)/2.0;
	if(debug_level>3)_DBG_<<"istart="<<istart<<" iend="<<iend<<" theta_min="<<seed.theta_min<<" theta_max="<<seed.theta_max<<endl;
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
	double tan_alpha_min = tan(theta_min - theta_err)/r0;
	double tan_alpha_max = tan(theta_max + theta_err)/r0;
	for(unsigned int i=0; i<seed.hits.size(); i++){
		DFDCTrkHit *trkhit = seed.hits[i];
		if(!trkhit->flags&VALID_HIT)continue;
		
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
		if(debug_level>3)_DBG_<<" -- phi_hit="<<trkhit->phi_hit<<" z_hit="<<z_hit<<endl;
		if(debug_level>3)_DBG_<<" -- zmin="<<zmin<<"  zmax="<<zmax<<endl;
		
		// Find index of bins corresponding to tmin and tmax
		unsigned int imin = (unsigned int)floor((zmin-hist_low_limit)/bin_width);
		unsigned int imax = (unsigned int)floor((zmax-hist_low_limit)/bin_width);
		
		// If entire range of this hit is outside of the histogram limit
		// then discard this hit.
		if(imax<=0 || imin>=Nbins)continue;
		
		// Clip limits of bin range to our histogram limits
		if(imin<0)imin=0;
		if(imin>=Nbins)imin=Nbins-1;
		if(imax<0)imax=0;
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
	if(hist[istart]==0.0)seed.valid=false;
	
	// Calculate z limits
	seed.z_min = hist_low_limit + bin_width*(0.5+(double)istart);
	seed.z_max = hist_low_limit + bin_width*(0.5+(double)iend);
	seed.z_vertex = (seed.z_max + seed.z_min)/2.0;
	if(debug_level>3)_DBG_<<"istart="<<istart<<" iend="<<iend<<" z_min="<<seed.z_min<<" z_max="<<seed.z_max<<" hits[istart]="<<hist[istart]<<endl;
}

#if 0

//------------------
// TraceSeed
//------------------
int DTrackCandidate_factory_FDC::TraceSeed(Dtrk_hit *hit)
{
	// Starting with "hit", look for nearest neighbors to
	// trace the seed out on one side as far as possible
	// until we find a hit that is either too far away, or
	// changes directions by too large an angle
	int N = 0;
	Dtrk_hit *current_hit = hit;
	Dtrk_hit *last_hit = NULL;
	do{
		N++;
		current_hit->flags |= Dtrk_hit::IN_SEED;
		seed.hits_in_seed.push_back(current_hit);
		if(seed.hits_in_seed.size() >= MAX_SEED_HITS)break;
		Dtrk_hit *hit = FindClosestXY(current_hit);
		if(!hit)break;
		if(hit->DistXY2(current_hit) > MAX_SEED_DIST2)break;
		if(last_hit){
			if(current_hit->CosPhiXY(hit,last_hit) > 0.0)break;
		}
				
		last_hit = current_hit;
		current_hit = hit;

	}while(current_hit);

	return N;
}

//------------------
// FindClosestXY
//------------------
Dtrk_hit* DTrackCandidate_factory_FDC::FindClosestXY(Dtrk_hit *hit)
{
	Dtrk_hit *closest_hit = NULL;
	unsigned int mask = Dtrk_hit::USED | Dtrk_hit::IN_SEED | Dtrk_hit::IGNORE;
	float d2_min = 1000.0;
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *a = trkhits[i];
		if(a->flags & mask)continue;
		float d2 = hit->DistXY2(a);
		if(d2<d2_min){
			closest_hit = a;
			d2_min = d2;
		}
	}
	
	return closest_hit;
}

//------------------
// FitSeed
//------------------
int DTrackCandidate_factory_FDC::FitSeed(void)
{
	// Do a quick circle fit to find the center of the seed
	for(unsigned int i=0; i<seed.hits_in_seed.size(); i++){
		Dtrk_hit *hit = seed.hits_in_seed[i];
		seed.seed_fit.AddHitXYZ(hit->X(), hit->Y(), hit->Z());
	}
	seed.seed_fit.FitCircle();
	x0 = seed.seed_fit.x0;
	y0 = seed.seed_fit.y0;
	r0 = sqrt(x0*x0 + y0*y0);
		
	// Set ON_CIRCLE flag for all hits close to the circle defined by x0,y0
	int N_in_seed_and_on_circle = 0;
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *hit = trkhits[i];
		hit->flags &= ~(Dtrk_hit::ON_CIRCLE | Dtrk_hit::ON_LINE | Dtrk_hit::ON_TRACK);
		if(hit->flags&Dtrk_hit::USED)continue;

		float dx = hit->X() - x0;
		float dy = hit->Y() - y0;
		float r = sqrt(dx*dx + dy*dy);
		
		if(DEBUG_HISTS){
			const DCDCTrackHit* cdchit = dynamic_cast<const DCDCTrackHit*>(eventLoop->FindByID(hit->hitid));
			if(cdchit){
				dist_to_seed_vs_cdclayer->Fill(cdchit->wire->ring, fabs(r0-r));
			}
		}
		
		if(fabs(r0-r) <= MAX_CIRCLE_DIST){
			hit->flags |= Dtrk_hit::ON_CIRCLE;
			seed.hits_on_seed_circle.push_back(hit);
			seed.circle_fit1.AddHitXYZ(hit->X(), hit->Y(), hit->Z());
			if(hit->flags & Dtrk_hit::IN_SEED)N_in_seed_and_on_circle++;
		}
	}

	// If the IN_SEED hits don't actually end up on the circle
	// (I don't quite understand it, but it happens) then we should
	// set the ignore flag for the first IN_SEED hit and try again.
	if(N_in_seed_and_on_circle<3 || seed.hits_on_seed_circle.size()<4){
		ChopSeed();
		return 0;
	}

	// Do a second circle fit using all of the on_seed_circle points
	seed.circle_fit1.FitCircle();
	x0 = seed.circle_fit1.x0;
	y0 = seed.circle_fit1.y0;
	r0 = sqrt(x0*x0 + y0*y0);

	// Now, make one more pass to flag the hits that are close to the circle fit
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *hit = trkhits[i];
		hit->flags &= ~(Dtrk_hit::ON_CIRCLE | Dtrk_hit::ON_LINE | Dtrk_hit::ON_TRACK);
		if(hit->flags&Dtrk_hit::USED)continue;

		float dx = hit->X() - x0;
		float dy = hit->Y() - y0;
		float r = sqrt(dx*dx + dy*dy);
			
		if(fabs(r0-r) <= MAX_CIRCLE_DIST){
			hit->flags |= Dtrk_hit::ON_CIRCLE;
			seed.hits_on_circle.push_back(hit);
			seed.circle_fit2.AddHitXYZ(hit->X(), hit->Y(), hit->Z());
		}
	}
	seed.circle_fit2.FitCircle();
	x0 = seed.circle_fit2.x0;
	y0 = seed.circle_fit2.y0;
	r0 = sqrt(x0*x0 + y0*y0);

	return 1;
}

//------------------
// FindLineHits
//------------------
int DTrackCandidate_factory_FDC::FindLineHits(void)
{
	// The hits_on_circle vector should now contain pointers to
	// only those hits on the circle. The hits_on_circle_with_z
	// vector should contain hits from CDC stereo wires that have
	// valid z coordinates and are on the circle. Add the FDC hits
	// from hits_on_circle to the hits from hits_on_circle_with_z
	// to get a list of all hits with valid z coordinate info
	// (i.e. space points).
	for(unsigned int i=0; i<seed.hits_on_circle.size(); i++){
		Dtrk_hit *hit = seed.hits_on_circle[i];
		if(hit->system == SYS_FDC)seed.hits_on_circle_with_z.push_back(hit);
	}
	
	// Order hits by z.
	sort(seed.hits_on_circle_with_z.begin(), seed.hits_on_circle_with_z.end(), TrkHitZSort_C);

	// Find the phi-z angle and the z-vertex
	int ok = FindPhiZAngle();
	if(ok)ok = FindZvertex();	
	if(!ok){
		ChopSeed();
		return 0;
	}

	// Find all on-circle hits consistent with phizangle and zvertex
	float cos_phizangle = cos(seed.phizangle);
	float sin_phizangle = sin(seed.phizangle);
	for(unsigned int i=0; i<seed.hits_on_circle_with_z.size(); i++){
		Dtrk_hit *hit = seed.hits_on_circle_with_z[i];
		float dphi = hit->phi_circle;
		float dz = hit->Z()-seed.z_vertex;
		float dr = sqrt(dphi*dphi + dz*dz);
		float sin_rel = (dphi*cos_phizangle - dz*sin_phizangle)/dr;
		float d = sin_rel*dr;
		if(fabs(d)<MAX_PHI_Z_DIST){
			// Flags hits as "ON_LINE" and push onto hits_on_line vector
			hit->flags |= Dtrk_hit::ON_LINE;
			seed.hits_on_line.push_back(hit);
		}
	}

	if(seed.hits_on_line.size()<MIN_PHI_Z_HITS){
		ChopSeed();
		return 0;
	}

	return 1;
}

//------------------
// FindPhiZAngle
//------------------
int DTrackCandidate_factory_FDC::FindPhiZAngle(void)
{
	// Fill the phi_circle field for all hits on the circle
	// centered at x0,y0
	Fill_phi_circle(seed.hits_on_circle_with_z, x0, y0);

	// Loop over hits on circle with z information and fill the phizangle histo.
	seed.phizangle_hist->Reset();
	for(unsigned int i=0; i<seed.hits_on_circle_with_z.size(); i++){
		Dtrk_hit *hit = seed.hits_on_circle_with_z[i];

		float theta1 = atan2f((hit->phi_circle - 0.0), hit->Z() - TARGET_Z_MIN);
		float theta2 = atan2f((hit->phi_circle - 0.0), hit->Z() - TARGET_Z_MAX);
		float theta_min = theta1<theta2 ? theta1:theta2;
		float theta_max = theta1<theta2 ? theta2:theta1;
		
		for(float t=theta_min; t<=theta_max; t+=phizangle_bin_size){
			seed.phizangle_hist->Fill(t);
		}
	}
	
	// Use maximum bin for phi-z angle
	int xbin, ybin, zbin;
	seed.phizangle_hist->GetMaximumBin(xbin,ybin,zbin);
	seed.phizangle = seed.phizangle_hist->GetXaxis()->GetBinCenter(xbin);
	
	// Bin content should be at least as musch as the minimum
	// number of hits required for a track.
	float height = seed.phizangle_hist->GetBinContent(xbin);
	if(height<3)return 0;

	// For large angle tracks, the "peak" in the phizangle histo
	// becomes a plateau. In this case, the maximum bin returned
	// is the left-most bin of the plateau. (At least this  is the
	// current behaviour of ROOT.) Look for a plateau and use the
	// center of it if it exists.
	int Nbins = seed.phizangle_hist->GetXaxis()->GetNbins();
	int xbin_right = xbin;
	for(int i=xbin+1; i<Nbins; i++){
		if(seed.phizangle_hist->GetBinContent(i) != height)break;
		xbin_right = i;
	}
	phizangle_min = seed.phizangle - phizangle_bin_size/2.0;
	phizangle_max = phizangle_min + phizangle_bin_size*(float)(xbin_right-xbin+1);
	seed.phizangle = (phizangle_min + phizangle_max)/2.0;

	return 1;
}

//------------------
// Fill_phi_circle
//------------------
void DTrackCandidate_factory_FDC::Fill_phi_circle(vector<Dtrk_hit*> &hits, float X, float Y)
{
	/// Fill in the phi_circle attribute for all Dtrk_hits in "hits" by
	/// calculating phi measured about the axis at x0,y0. Hits with either
	/// the IGNORE or USED flags set will be ignored.
	float r0 = sqrt(X*X + Y*Y);
	float x_last = -X;
	float y_last = -Y;
	float phi_last = 0.0;
	for(unsigned int i=0; i<hits.size(); i++){
		Dtrk_hit *hit = hits[i];

		float dx = hit->X() - X;
		float dy = hit->Y() - Y;
		float dphi = atan2f(dx*y_last - dy*x_last, dx*x_last + dy*y_last);
		float phi = phi_last +dphi;
		hit->phi_circle = phi*r0;
		x_last = dx;
		y_last = dy;
		phi_last = phi;
	}
}

//------------------
// FindZvertex
//------------------
int DTrackCandidate_factory_FDC::FindZvertex(void)
{
	// Reset z_vertex histogram
	seed.zvertex_hist->Reset();

	// tan(x) has poles at x=+/- pi/2 What we really want
	// is cot(x) but that isn't always available. To avoid
	// the pole condition, calculate it using cos(x)/sin(x).
	float cot_phizangle_min = cos(phizangle_min)/sin(phizangle_min);
	float cot_phizangle_max = cos(phizangle_max)/sin(phizangle_max);

	// Loop over hits with z information
	for(unsigned int i=0; i<seed.hits_on_circle_with_z.size(); i++){
		Dtrk_hit *hit = seed.hits_on_circle_with_z[i];

		// Find intersections with Z-axis for lines with min and
		// max phizangle going through point (a->z, a->phi_circle).
		float z1 = hit->Z() - hit->phi_circle*cot_phizangle_min;
		float z2 = hit->Z() - hit->phi_circle*cot_phizangle_max;
		float z_min = z1<z2 ? z1:z2;
		float z_max = z1<z2 ? z2:z1;

		if(z_min<TARGET_Z_MIN) z_min = TARGET_Z_MIN;
		if(z_max>TARGET_Z_MAX) z_max = TARGET_Z_MAX;
		
		for(float z=z_min; z<=z_max; z+=z_vertex_bin_size){
			seed.zvertex_hist->Fill(z);
		}
	}

	// Use maximum bin for z vertex
	int xbin, ybin, zbin;
	seed.zvertex_hist->GetMaximumBin(xbin,ybin,zbin);
	seed.z_vertex = seed.zvertex_hist->GetXaxis()->GetBinCenter(xbin);
	
	// Bin content should be at least as musch as the minimum
	// number of hits required for a track.
	float height = seed.zvertex_hist->GetBinContent(xbin);
	if(height<3)return 0;

	// Similar to phizangle method above, find range of plateau
	int Nbins = seed.zvertex_hist->GetXaxis()->GetNbins();
	int xbin_right = xbin;
	for(int i=xbin+1; i<Nbins; i++){
		if(seed.zvertex_hist->GetBinContent(i) != height)break;
		xbin_right = i;
	}
	float z_vertex_min = seed.z_vertex - z_vertex_bin_size/2.0;
	float z_vertex_max = z_vertex_min + z_vertex_bin_size*(float)(xbin_right-xbin+1);
	seed.z_vertex = (z_vertex_min + z_vertex_max)/2.0;

	return 1;
}


//------------------
// FitTrack
//------------------
int DTrackCandidate_factory_FDC::FitTrack(void)
{
	// At this point, the values in x0, y0 come from a circle fit
	// (linear regression) involving CDC axial and FDC hits. It's
	// possible some of the FDC hits do not fall on the phi-z
	// line and so, should be excluded from the circle fit.
	//
	// Also, the values of z_vertex and phizangle come from histograms
	// and so are restricted by the bin widths. These are filled by
	// CDC stereo and FDC hits. Because there are often only a few points
	// to fit a line to in the phi-z plane, a linear regression will
	// often give a result that does not point back to the target.
	// We can however, still get a better value for theta (derived
	// from phizangle) by doing a fit that restricts the z_vertex
	// to being the current value.
	//
	// Since both of these "fits" use different (although overlapping)
	// data sets, we cannot really do them simultaneously with DQuickFit.
	// Thus, the two fits below.
	
	// For now we just use the last circle fit and phizangle, z_vertex
	// to decide the parameters. 
	DQuickFit &fit = seed.circle_fit2;

	double phi = fit.phi;
	double p_trans = fit.p_trans;

	// Fit to line to get z_vertex and theta
	seed.line_fit = seed.circle_fit2;
	seed.line_fit.Clear();
	for(unsigned int i=0; i<seed.hits_on_line.size(); i++){
		Dtrk_hit *a = seed.hits_on_line[i];
		seed.line_fit.AddHitXYZ(a->X(), a->Y(), a->Z());
	}
	seed.line_fit.FitLine_FixedZvertex(seed.z_vertex); // circle parameters from fit above are still valid
	
	// The following is from a fit of ratio of thrown to reconstructed
	// transverse momentum vs. theta
	double par[] = {0.984463, 0.150759, -0.414933, 0.257472, -0.055801};
	double theta = seed.line_fit.theta;
	double ff = par[0]+theta*(par[1]+theta*(par[2]+theta*(par[3]+theta*par[4])));
	double sin_theta = sin(theta);

	// Create a new DTrackCandidate object
	DTrackCandidate *can = new DTrackCandidate;

	can->x0 = fit.x0;
	can->y0 = fit.y0;
	can->p_trans = p_trans*ff;
	can->phi = phi;
	can->z_vertex = seed.z_vertex;
	can->dzdphi = r0/theta;
	can->p = can->p_trans/sin_theta;
	can->q = seed.line_fit.q;
	can->theta = theta;

	// Fill in DKinematic Data protion of this
	can->setMass(0.0);
	can->setMomentum(DVector3(can->p_trans*cos(can->phi), can->p_trans*sin(can->phi), can->p*cos(can->theta)));
	can->setPosition(DVector3(0.0, 0.0, can->z_vertex));
	can->setCharge(can->q);

	// Do one last pass over the hits, marking all 
	// that are consistent with these parameters as used.
	MarkTrackHits(can, &seed.line_fit);
	
	_data.push_back(can);

	return 1;
}

//------------------
// MarkTrackHits
//------------------
int DTrackCandidate_factory_FDC::MarkTrackHits(DTrackCandidate *trackcandidate, DQuickFit *fit)
{
	// Mark all hits_on_line and hits_on_circle as used and simultaneously
	// fill the hitid member of the DrackCandidate object
	int N_in_seed_and_on_track = 0;
	for(unsigned int i=0; i<seed.hits_on_circle.size(); i++){
		Dtrk_hit *hit = seed.hits_on_circle[i];
		hit->flags |= Dtrk_hit::ON_TRACK | Dtrk_hit::USED;
		if(hit->flags & Dtrk_hit::IN_SEED)N_in_seed_and_on_track++;

		// Add hit index to track in factory data
		trackcandidate->hitid.push_back(hit->hitid);
		seed.hits_on_track.push_back(hit);
	}
	
	// It's possible that none of the seed hits are included in the
	// final track, when that happens, we have to ignore the first
	// seed hit on the next loop or else we'll keep getting the same
	// seed and the same result.
	if(seed.hits_on_track.size()<3 || N_in_seed_and_on_track<1){
		ChopSeed();
		return 0;
	}

	return 1;
}

#endif

//------------------
// toString
//------------------
const string DTrackCandidate_factory_FDC::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("      id: Nhits: q:     p:       theta:   phi:  p_trans:   x:     y:     z:    dz/dphi:");

	for(unsigned int i=0; i<_data.size(); i++){
		DTrackCandidate *trackcandidate = _data[i];
		printnewrow();
		
		printcol("%lx",    trackcandidate->id);
#if 0
		printcol("%d",    trackcandidate->hitid.size());
		printcol("%+d", (int)trackcandidate->q);
		printcol("%3.3f", trackcandidate->p);
		printcol("%1.3f", trackcandidate->theta);
		printcol("%1.3f", trackcandidate->phi);
		printcol("%3.2f", trackcandidate->p_trans);
		printcol("%2.2f", trackcandidate->x0);
		printcol("%2.2f", trackcandidate->y0);
		printcol("%2.2f", trackcandidate->z_vertex);
		printcol("%1.3f", trackcandidate->dzdphi);
#endif
		printrow();
	}
	
	return _table;
}


// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ----------------------- DEBUGGING ROUTINES ---------------------
// ----------------------------------------------------------------
// ----------------------------------------------------------------

#if 0
//------------------
// DumpHits
//------------------
void DTrackCandidate_factory_FDC::DumpHits(int current_seed_number, string stage)
{
	// Count statistics for flags
	int Nhits = 0;
	int in_seed=0, used=0, on_circle=0, ignore=0 ,on_line=0, on_track=0;
	for(unsigned int i=0; i<trkhits.size(); i++){
		const Dtrk_hit *a = trkhits[i];
		Nhits++;
		if(a->flags&Dtrk_hit::IN_SEED)in_seed++;
		if(a->flags&Dtrk_hit::USED)used++;
		if(a->flags&Dtrk_hit::ON_CIRCLE)on_circle++;
		if(a->flags&Dtrk_hit::IGNORE)ignore++;
		if(a->flags&Dtrk_hit::ON_LINE)on_line++;
		if(a->flags&Dtrk_hit::ON_TRACK)on_track++;
	}

	// Print message
	cout<<endl;
	cout<<"----------- DUMPING TRACK CANDIDATE HITS ------------"<<endl;
	cout<<"RUN / EVENT: "<<runnumber<<" / "<<eventnumber<<endl;
	cout<<"SEED NUMBER: "<<current_seed_number<<endl;
	cout<<"      STAGE: "<<stage<<endl;
	cout<<"    IN_SEED: "<<in_seed<<endl;
	cout<<"  ON_CIRCLE: "<<on_circle<<endl;
	cout<<"    ON_LINE: "<<on_line<<endl;
	cout<<"   ON_TRACK: "<<on_track<<endl;
	cout<<"       USED: "<<used<<endl;
	cout<<"     IGNORE: "<<ignore<<endl;
	
	// Print more verbose info
	if(debug_level<5)return;
	cout<<"hits:"<<endl;
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *t = trkhits[i];
		
		cout<<i<<" flags:"<<t->flags<<" = ";
		if(t->flags == 0x0)cout<<" NONE";
		if(t->flags & Dtrk_hit::USED)cout<<" USED";
		if(t->flags & Dtrk_hit::IN_SEED)cout<<" IN_SEED";
		if(t->flags & Dtrk_hit::ON_CIRCLE)cout<<" ON_CIRCLE";
		if(t->flags & Dtrk_hit::ON_LINE)cout<<" ON_LINE";
		if(t->flags & Dtrk_hit::ON_TRACK)cout<<" ON_TRACK";
		if(t->flags & Dtrk_hit::IGNORE)cout<<" IGNORE";
		cout<<endl;
		cout<<i<<" phi_circle: "<<t->phi_circle<<endl;
		cout<<i<<" x y z: "<<t->X()<<" "<<t->Y()<<" "<<t->Z()<<endl;
		cout<<i<<" r phi: "<<t->Mag()<<" "<<t->Phi()<<endl;
		cout<<i<<" system:"<<t->system;
		if(t->system == 0x0)cout<<" NONE";
		if(t->system & SYS_CDC)cout<<" CDC";
		if(t->system & SYS_FDC)cout<<" FDC";
		if(t->system & SYS_BCAL)cout<<" BCAL";
		if(t->system & SYS_FCAL)cout<<" FCAL";
		if(t->system & SYS_TOF)cout<<" TOF";
		if(t->system & SYS_UPV)cout<<" UPV";
		if(t->system & SYS_CHERENKOV)cout<<" CHERENKOV";
		if(t->system & SYS_TAGGER)cout<<" TAGGER";
		cout<<endl;
		cout<<endl;		
	}
}

#endif
