// $Id$
//
//    File: DTrackCandidate_factory.cc
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>

#include <TH1.h>

#include "DTrackCandidate_factory.h"
#include "JANA/JApplication.h"
#include "DTrack.h"
#include "DQuickFit.h"
#include "JANA/JGeometry.h"
#include "DMagneticFieldMap.h"
#include "DVector2.h"


bool TrkHitSort_C(Dtrk_hit* const &hit1, Dtrk_hit* const &hit2) {
	return hit1->Mag() > hit2->Mag();
}

bool TrkHitZSort_C(Dtrk_hit* const &hit1, Dtrk_hit* const &hit2) {
	return hit1->Z() < hit2->Z();
}


//------------------
// DTrackCandidate_factory
//------------------
DTrackCandidate_factory::DTrackCandidate_factory()
{
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
	
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_DIST",		MAX_SEED_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_SEED_HITS",		MAX_SEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_SEED_HITS",		MIN_SEED_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_CIRCLE_DIST",		MAX_CIRCLE_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_PHI_Z_DIST",		MAX_PHI_Z_DIST);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_PHI_Z_HITS",		MIN_PHI_Z_HITS);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_DEBUG_BUFFERS",	MAX_DEBUG_BUFFERS);
	gPARMS->SetDefaultParameter("TRKFIND:TARGET_Z_MIN",			TARGET_Z_MIN);
	gPARMS->SetDefaultParameter("TRKFIND:TARGET_Z_MAX",			TARGET_Z_MAX);
	gPARMS->SetDefaultParameter("TRKFIND:TRACKHIT_SOURCE",		TRACKHIT_SOURCE);
	gPARMS->SetDefaultParameter("TRKFIND:XY_NOISE_CUT",			XY_NOISE_CUT);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_HIT_Z",				MIN_HIT_Z);
	gPARMS->SetDefaultParameter("TRKFIND:MAX_HIT_Z",				MAX_HIT_Z);
	gPARMS->SetDefaultParameter("TRKFIND:EXCLUDE_STEREO",		EXCLUDE_STEREO);
	gPARMS->SetDefaultParameter("TRKFIND:MIN_CANDIDATE_HITS",	MIN_CANDIDATE_HITS);

	MAX_SEED_DIST2 = MAX_SEED_DIST*MAX_SEED_DIST;
	XY_NOISE_CUT2 = XY_NOISE_CUT*XY_NOISE_CUT;
	
}

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory::init(void)
{
	dgeom = NULL;
	bfield = NULL;
	
	char suffix[32];
	sprintf(suffix,"_%08x", (unsigned int)pthread_self());
	
	// Since there is a subclass of this one (DTrackCandidate_factory_THROWN)
	// that does not need these histograms, we check if we are a true 
	// DTrackCandidate_factory object or not before creating them. If we
	// don't then the DTrackCandidate_factory_THROWN class will recreate
	// them leading to warning messages from ROOT.
	if(typeid(this) == typeid(DTrackCandidate_factory*)){
		char name[64];
		eventLoop->GetJApplication()->Lock();

		sprintf(name,"phi_z_angle%s",suffix);
		seed.phizangle_hist = new TH1F(name,"phi_z_angle", 1000, -M_PI, M_PI);
		phizangle_bin_size = seed.phizangle_hist->GetBinCenter(2) - seed.phizangle_hist->GetBinCenter(1);
		
		sprintf(name,"z_vertex%s",suffix);
		seed.zvertex_hist = new TH1F(name,"z_vertex", 140, TARGET_Z_MIN, TARGET_Z_MAX);
		z_vertex_bin_size = seed.zvertex_hist->GetBinCenter(2) - seed.zvertex_hist->GetBinCenter(1);

		eventLoop->GetJApplication()->Unlock();
	}
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory::brun(JEventLoop *loop, int runnumber)
{
	this->runnumber = runnumber;
	dgeom = loop->GetJApplication()->GetJGeometry(runnumber);
	bfield = NULL;
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory::fini(void)
{
	// seems errors are sometimes caused by deleting these ???
	eventLoop->GetJApplication()->Lock();
	if(seed.phizangle_hist)delete seed.phizangle_hist;
	if(seed.zvertex_hist)delete seed.zvertex_hist;
	eventLoop->GetJApplication()->Unlock();
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory::evnt(JEventLoop *loop, int eventnumber)
{
	this->eventnumber = eventnumber;

	// Clear previous event from internal buffers
	ClearEvent();

	// Get the hits into the trkhits vector
	GetTrkHits(loop);

	// Loop as long as we keep finding tracks
	for(int i=0;i<100;i++){
		// Clear vectors of various flavors of hits
		seed.Clear();
	
		// Find a seed (group of hits that appear to be a clean track segment)
		if(!FindSeed()){
			if(debug_level>3)DumpHits(i, "FindSeed");
			break;
		}
		
		// Fit the seed hits to a circle and flag all unused hits
		// which fall close to the circle. If the fit fails, the
		// first IN_SEED hit automatically has it's IGNORE flag set
		// so we can just jump to the next iteration of the loop.
		if(!FitSeed()){
			if(debug_level>3)DumpHits(i, "FitSeed");
			dbg_seeds.push_back(seed);
			continue;
		}

		// The CDC hits do not have valid z-info at this point. Use
		// the stereo layers to generate it and at the same time, pull
		// in the stereo hits.
		if(!FindCDCStereoZvals()){
			if(debug_level>3)DumpHits(i,"FindCDCStereoZvals");
			dbg_seeds.push_back(seed);
			continue;
		}
		
		// Using results of X/Y fit, find phi-z angle and z_vertex
		// Make a list of hits consistent with the values found.
		if(!FindLineHits()){
			if(debug_level>3)DumpHits(i,"FindLineHits");
			dbg_seeds.push_back(seed);
			continue;
		}
		
		// Fit the hits in the list created by above. Use result
		// to make a new DTrackCandidate
		FitTrack();
		if(debug_level>3)DumpHits(i,"FitTrack");
		dbg_seeds.push_back(seed);
	}

#if 0
	// Filter out "bad" track candidates
	for(unsigned int i=0; i<_data.size(); i++){
		DTrackCandidate *trackcandidate = _data[i];

		if(trackcandidate->hitid.size()<MIN_CANDIDATE_HITS){
			delete trackcandidate;
			_data.erase(_data.begin()+i);
			if(dbg_track_fit.size()<MAX_DEBUG_BUFFERS){
				delete dbg_track_fit[i];
				dbg_track_fit.erase(dbg_track_fit.begin()+i);
				dbg_hol.erase(dbg_hol.begin()+i);
				dbg_phizangle.erase(dbg_phizangle.begin()+i);
				dbg_z_vertex.erase(dbg_z_vertex.begin()+i);
				dbg_seed_index.erase(dbg_seed_index.begin()+i);
			}
			i--;
		}
	}
#endif

	return NOERROR;
}

//------------------
// ClearEvent
//------------------
void DTrackCandidate_factory::ClearEvent(void)
{
	// Clear TrkHits from previous event
	for(unsigned int i=0; i<trkhits.size(); i++)delete trkhits[i];
	trkhits.clear();
	for(unsigned int i=0; i<trkhits_stereo.size(); i++)delete trkhits_stereo[i];
	trkhits_stereo.clear();
	for(unsigned int i=0; i<trkhits_extra.size(); i++)delete trkhits_extra[i];
	trkhits_extra.clear();

	// Clear debugging info (if any) from previous event
	dbg_seeds.clear();

}

//------------------
// GetTrkHits
//------------------
void DTrackCandidate_factory::GetTrkHits(JEventLoop *loop)
{
	// Copy Hits into trkhits vector. The objects returned as
	// DTrackHit should really be Dtrk_hit types. Dynamically
	// cast them back
	cdctrackhits.clear();
	fdcpseudos.clear();
	loop->Get(cdctrackhits);
	loop->Get(fdcpseudos);

	// Create Dtrk_hit objects from CDC hits
	for(unsigned int i=0; i<cdctrackhits.size(); i++){
		const DCDCTrackHit *cdchit = cdctrackhits[i];

		Dtrk_hit *hit = new Dtrk_hit(cdchit->wire->origin);
		hit->flags = 0x0;
		hit->hitid = cdchit->id;
		hit->wire = cdchit->wire;
		hit->system = SYS_CDC;

		// Need to Fill covariance matrix for hit
		
		if(cdchit->wire->stereo==0.0){
			// Axial wire
			trkhits.push_back(hit);
		}else{
			// Stereo wire
			trkhits_stereo.push_back(hit);
		}
	}

	// Create Dtrk_hit objects from FDC hits
	// Note that these have already undergone some level of reconstruction
	for(unsigned int i=0; i<fdcpseudos.size(); i++){
		const DFDCPseudo *fdcpseudo = fdcpseudos[i];
		
		DVector3 pos = fdcpseudo->wire->origin + fdcpseudo->s*fdcpseudo->wire->udir;
		Dtrk_hit *hit = new Dtrk_hit(pos);
		hit->flags = 0x0;
		hit->hitid = fdcpseudo->id;
		hit->wire = fdcpseudo->wire;
		hit->system = SYS_FDC;
		
		// Need to Fill covariance matrix for hit
		
		trkhits.push_back(hit);
	}

	
	// Sort hits by r in X/Y plane
	trkhits_r_sorted = trkhits;
	sort(trkhits_r_sorted.begin(), trkhits_r_sorted.end(), TrkHitSort_C);

	// Order the track hits by z.
	sort(trkhits.begin(), trkhits.end(), TrkHitZSort_C);

	// Flag all "lone" hits to be ignored
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *hit = trkhits[i];
		
		// temporarily flag this hit to be ignored so FindClosestXY will work
		hit->flags |= Dtrk_hit::IGNORE;
		Dtrk_hit *b = FindClosestXY(hit);
		if(b){
			if(hit->DistXY2(b) < XY_NOISE_CUT2)hit->flags &= ~Dtrk_hit::IGNORE;
		}
	}
}

//------------------
// FindSeed
//------------------
int DTrackCandidate_factory::FindSeed(void)
{

	// Loop over all un-used, non-ignored hits and try
	// tracing a seed from each. Once a seed with MIN_SEED_HITS
	// or more hits is found, return a 1 so a track can
	// be searched for.
	for(unsigned int i=0; i<trkhits_r_sorted.size(); i++){
		Dtrk_hit *hit = trkhits_r_sorted[i];
		if(!(hit->flags & (Dtrk_hit::USED | Dtrk_hit::IGNORE))){
		
			// Clear IN_SEED bit flag all hits
			for(unsigned int j=0; j<trkhits.size(); j++){
				trkhits[j]->flags &= ~(Dtrk_hit::IN_SEED);
			}

			// Trace the seed out by looking for nearest neighbors.
			int N = TraceSeed(hit);
			if(N>=(int)MIN_SEED_HITS)return 1;
			ChopSeed();
		}
	}
	
	// Hmmm... looks like no seeds could be found. Return 0
	
	return 0;
}

//------------------
// TraceSeed
//------------------
int DTrackCandidate_factory::TraceSeed(Dtrk_hit *hit)
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
Dtrk_hit* DTrackCandidate_factory::FindClosestXY(Dtrk_hit *hit)
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
int DTrackCandidate_factory::FitSeed(void)
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
// FindCDCStereoZvals
//------------------
int DTrackCandidate_factory::FindCDCStereoZvals(void)
{
	/// Find the z-coordinates for the CDC hits coming from the stereo layers
	/// and add them to the list of hits.
	
	// To find the z coordinate, we look at the 2D projection of the
	// stereo wire and find the intersection point of that with the
	// circle found in FitSeed().

	// Loop over stereo wires to find z-values for each
	for(unsigned int i=0; i<trkhits_stereo.size(); i++){
		Dtrk_hit *hit = trkhits_stereo[i];
		const DCoordinateSystem *wire = hit->wire;
		DVector2 r1(wire->origin.X(), wire->origin.Y());
		DVector2 r2(wire->udir.X(), wire->udir.Y());
		DVector2 R(x0, y0);
		r2*=1.0/r2.Mod();
		double a = 1.0; // r2.Mod2()
		double b = 2.0*r2*(r1-R);
		double c = r1.Mod2()-2.0*r1*R;
		double A = b*b - 4.0*a*c;

		// Check that this wire intersects this circle
		if(A<0.0)continue; // line along wire does not intersect circle, ever.
		if(a==0.0){
			_DBG_<<"wire in CDC stereo hit list is not stereo!"<<endl;
			continue; // this must not be a stereo wire!
		}

		// Calculate intersection points
		double B = sqrt(A);
		double alpha1 = (-b - B)/(2.0*a);
		double alpha2 = (-b + B)/(2.0*a);
		
		// At this point we must decide which value of alpha to use. The
		// proper way would likely involve either trying both possibilities
		// and taking the one that gave the better chi-sq for a line fit
		// of phi vs. z or looking at the surrounding axial layers
		// and using the value which puts the hit closest to those.
		// For now, we just use the value closest to zero (i.e. closest to
		// the center of the wire).
		double alpha = fabs(alpha1)<fabs(alpha2) ? alpha1:alpha2;
		
		// Now we must convert the alpha value into a z-value. To do this,
		// we use the known theta angle of the wire. The distance along the
		// wire from the center in 3D is given by:
		//
		//   s = alpha/sin(theta_wire)
		//
		// The z coordinate of the hit is given by:
		//
		//    z = z1 + s*z2
		//
		// where z1 is the z coordinate of the center of the wire and z2
		// is the z component of the direction of the wire. i.e.
		//
		//    z2 = cos(theta_wire)
		//
		// This means  sin(theta_wire) = sqrt(1 - (z2)^2)
		double z2 = wire->udir.Z();
		double s = alpha/sqrt(1.0-z2*z2);
		if(fabs(s) > wire->L)continue; // if wire doesn't cross circle, skip hit
		
		// Create a new hit object so we can set x,y,z based on s
		Dtrk_hit *hitz = new Dtrk_hit(wire->origin + s*wire->udir, *hit);
		trkhits_extra.push_back(hitz); // keep track of object for later deletion
		seed.hits_on_circle.push_back(hitz);
		seed.hits_on_circle_with_z.push_back(hitz);
	}


	return 1;
}

//------------------
// FindLineHits
//------------------
int DTrackCandidate_factory::FindLineHits(void)
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
int DTrackCandidate_factory::FindPhiZAngle(void)
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
void DTrackCandidate_factory::Fill_phi_circle(vector<Dtrk_hit*> &hits, float X, float Y)
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
int DTrackCandidate_factory::FindZvertex(void)
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
int DTrackCandidate_factory::FitTrack(void)
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
	DTrackCandidate *trackcandidate = new DTrackCandidate;

	trackcandidate->p_trans = p_trans*ff;
	trackcandidate->phi = phi;
	trackcandidate->z_vertex = seed.z_vertex;
	trackcandidate->dzdphi = r0/theta;
	trackcandidate->p = trackcandidate->p_trans/sin_theta;
	trackcandidate->q = seed.line_fit.q;
	trackcandidate->theta = theta;

	// Do one last pass over the hits, marking all 
	// that are consistent with these parameters as used.
	MarkTrackHits(trackcandidate, &seed.line_fit);
	
	_data.push_back(trackcandidate);

	return 1;
}

//------------------
// MarkTrackHits
//------------------
int DTrackCandidate_factory::MarkTrackHits(DTrackCandidate *trackcandidate, DQuickFit *fit)
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


//------------------
// toString
//------------------
const string DTrackCandidate_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("      id: Nhits: q:     p:       theta:   phi:  p_trans:   x:     y:     z:    dz/dphi:");

	for(unsigned int i=0; i<_data.size(); i++){
		DTrackCandidate *trackcandidate = _data[i];
		printnewrow();
		
		printcol("%lx",    trackcandidate->id);
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

		printrow();
	}
	
	return _table;
}


// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ----------------------- DEBUGGING ROUTINES ---------------------
// ----------------------------------------------------------------
// ----------------------------------------------------------------

//------------------
// DumpHits
//------------------
void DTrackCandidate_factory::DumpHits(int current_seed_number, string stage)
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


