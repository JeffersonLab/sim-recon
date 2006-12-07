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


class TrkHitSort{
	public:
		bool operator()(Dtrk_hit* const &hit1, Dtrk_hit* const &hit2) const {
			return hit1->r > hit2->r;
		}
};
class TrkHitZSort{
	public:
		bool operator()(Dtrk_hit* const &hit1, Dtrk_hit* const &hit2) const {
			return hit1->z < hit2->z;
		}
};

bool TrkHitSort_C(Dtrk_hit* const &hit1, Dtrk_hit* const &hit2) {
	return hit1->r > hit2->r;
}

bool TrkHitZSort_C(Dtrk_hit* const &hit1, Dtrk_hit* const &hit2) {
	return hit1->z < hit2->z;
}


//------------------
// DTrackCandidate_factory
//------------------
DTrackCandidate_factory::DTrackCandidate_factory()
{
	// Set defaults
	MAX_SEED_DIST = 5.0;
	MAX_SEED_HITS = 10;
	MAX_CIRCLE_DIST = 2.0;
	MAX_PHI_Z_DIST = 10.0;
	MAX_DEBUG_BUFFERS = 0;
	TARGET_Z_MIN = 50.0;
	TARGET_Z_MAX = 80.0;
	TRACKHIT_SOURCE = "MC";
	XY_NOISE_CUT = 2.0;
	MIN_HIT_Z = -100.0;
	MAX_HIT_Z = +360.0;
	EXCLUDE_STEREO = true;
	
	gPARMS->SetDefaultParameter("TRK:MAX_SEED_DIST",		MAX_SEED_DIST);
	gPARMS->SetDefaultParameter("TRK:MAX_SEED_HITS",		MAX_SEED_HITS);
	gPARMS->SetDefaultParameter("TRK:MAX_CIRCLE_DIST",	MAX_CIRCLE_DIST);
	gPARMS->SetDefaultParameter("TRK:MAX_PHI_Z_DIST",	MAX_PHI_Z_DIST);
	gPARMS->SetDefaultParameter("TRK:MAX_DEBUG_BUFFERS",MAX_DEBUG_BUFFERS);
	gPARMS->SetDefaultParameter("TRK:TARGET_Z_MIN",		TARGET_Z_MIN);
	gPARMS->SetDefaultParameter("TRK:TARGET_Z_MAX",		TARGET_Z_MAX);
	gPARMS->SetDefaultParameter("TRK:TRACKHIT_SOURCE",	TRACKHIT_SOURCE);
	gPARMS->SetDefaultParameter("TRK:XY_NOISE_CUT",		XY_NOISE_CUT);
	gPARMS->SetDefaultParameter("TRK:MIN_HIT_Z",			MIN_HIT_Z);
	gPARMS->SetDefaultParameter("TRK:MAX_HIT_Z",			MAX_HIT_Z);
	gPARMS->SetDefaultParameter("TRK:EXCLUDE_STEREO",	EXCLUDE_STEREO);

	MAX_SEED_DIST2 = MAX_SEED_DIST*MAX_SEED_DIST;
	XY_NOISE_CUT2 = XY_NOISE_CUT*XY_NOISE_CUT;
	
	dgeom = NULL;
	bfield = NULL;
	
	char suffix[32];
	sprintf(suffix,"_%08x", (unsigned int)pthread_self());
	
	char title[64];
	sprintf(title,"phi_z_angle%s",suffix);
	phizangle_hist = new TH1F(title,"phi_z_angle", 1000, -M_PI, M_PI);
	phizangle_bin_size = phizangle_hist->GetBinCenter(2)-phizangle_hist->GetBinCenter(1);
	sprintf(title,"phi_relative%s",suffix);
	phi_relative = new TH1F(title,"phi_relative", 1000, -M_PI, M_PI);
	sprintf(title,"z_vertex%s",suffix);
	zvertex_hist = new TH1F(title,"z_vertex", 140, TARGET_Z_MIN, TARGET_Z_MAX);
	z_vertex_bin_size = zvertex_hist->GetBinCenter(2)-zvertex_hist->GetBinCenter(1);
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory::brun(JEventLoop *loop, int runnumber)
{
	this->runnumber = runnumber;
	dgeom = loop->GetJApplication()->GetJGeometry(runnumber);
	bfield = NULL;
	//bfield = dgeom->GetDMagneticFieldMap();
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory::fini(void)
{
	// seems errors are sometimes caused by deleting these ???
	//delete phizangle_hist;
	//delete phi_relative;
	//delete zvertex_hist;
	
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
			continue;
		}
		
		// Using results of X/Y fit, find phi-z angle and z_vertex
		// Make a list of hits consistent with the values found.
		if(!FindLineHits()){
			if(debug_level>3)DumpHits(i,"FindLineHits");
			continue;
		}
		
		// Fit the hits in the list created by above. Use result
		// to make a new DTrackCandidate
		FitTrack();
		if(debug_level>3)DumpHits(i,"FitTrack");
	}
	
	// Filter out "bad" track candidates
	for(unsigned int i=0; i<_data.size(); i++){
		DTrackCandidate *trackcandidate = _data[i];

		if(trackcandidate->hitid.size()<10){
			delete trackcandidate;
			_data.erase(_data.begin()+i);
			if(dbg_track_fit.size()<MAX_DEBUG_BUFFERS){
				delete dbg_track_fit[i];
				dbg_track_fit.erase(dbg_track_fit.begin()+i);
			}
			i--;
		}
	}

	return NOERROR;
}

//------------------
// ClearEvent
//------------------
void DTrackCandidate_factory::ClearEvent(void)
{
	// Clear TrkHits from previous event
	trkhits.clear();

	// Clear debugging buffers
	if(MAX_DEBUG_BUFFERS){
		dbg_in_seed.clear();
		dbg_hoc.clear();
		dbg_hol.clear();
		dbg_hot.clear();
		for(unsigned int i=0; i<dbg_seed_fit.size(); i++)delete dbg_seed_fit[i];
		dbg_seed_fit.clear();
		for(unsigned int i=0; i<dbg_track_fit.size(); i++)delete dbg_track_fit[i];
		dbg_track_fit.clear();
		dbg_seed_index.clear();
		for(unsigned int i=0; i<dbg_phiz_hist.size(); i++)delete dbg_phiz_hist[i];
		dbg_phiz_hist.clear();
		dbg_phiz_hist_seed.clear();
		for(unsigned int i=0; i<dbg_zvertex_hist.size(); i++)delete dbg_zvertex_hist[i];
		dbg_zvertex_hist.clear();
		dbg_zvertex_hist_seed.clear();
		dbg_phizangle.clear();
		dbg_z_vertex.clear();
	}
}

//------------------
// GetTrkHits
//------------------
void DTrackCandidate_factory::GetTrkHits(JEventLoop *loop)
{
	// Copy Hits into trkhits vector. The objects returned as
	// DTrackHit should really be Dtrk_hit types. Dynamically
	// cast them back
	vector<const DTrackHit*> trackhits;
	loop->Get(trackhits, TRACKHIT_SOURCE.c_str());
	for(unsigned int i=0; i<trackhits.size(); i++){
		Dtrk_hit *hit = (Dtrk_hit*)dynamic_cast<const Dtrk_hit*>(trackhits[i]);
		if(hit){
			if(hit->system & (SYS_CDC|SYS_FDC)){
				if(hit->z>=MIN_HIT_Z && hit->z<=MAX_HIT_Z){
					// Don't include hits from stereo layers in CDC
					if(EXCLUDE_STEREO){
						if(hit->system==SYS_CDC){
							if(hit->r>22.5 && hit->r<31.5)continue;
							if(hit->r>40.5 && hit->r<48.5)continue;
						}
					}
					trkhits.push_back(hit);
				}
			}
		}
	}
	
	// Sort hits by r in X/Y plane
	trkhits_r_sorted = trkhits;
	sort(trkhits_r_sorted.begin(), trkhits_r_sorted.end(), TrkHitSort_C);

	// Order the track hits by z.
	sort(trkhits.begin(), trkhits.end(), TrkHitZSort_C);

	// Flag all "lone" hits to be ignored
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *a = trkhits[i];
		
		// temporarily flag this hit to be ignored so FindClosestXY will work
		a->flags |= Dtrk_hit::IGNORE;
		Dtrk_hit *b = FindClosestXY(a);
		if(b){
			if(a->DistXY2(b) < XY_NOISE_CUT2)a->flags &= ~Dtrk_hit::IGNORE;
		}
	}
}

//------------------
// FindSeed
//------------------
int DTrackCandidate_factory::FindSeed(void)
{

	// Loop over all un-used, non-ignored hits and try
	// tracing a seed from each. Once a seed with 4 or
	// more hits is found, return a 1 so a track can
	// be searched for.
	for(unsigned int i=0; i<trkhits_r_sorted.size(); i++){
		Dtrk_hit *a = trkhits_r_sorted[i];
		if(!(a->flags & (Dtrk_hit::USED | Dtrk_hit::IGNORE))){
		
			// Clear IN_SEED bit flag all hits
			for(unsigned int j=0; j<trkhits.size(); j++){
				trkhits[j]->flags &= ~(Dtrk_hit::IN_SEED);
			}

			// We call TraceSeed twice here on purpose! The first
			// call traces out in one direction and the second
			// will trace in the other.
			// NOTE: This will have a problem if hit "a"
			// is at a kink.
			hits_in_seed.clear();
			int N = TraceSeed(a);
			N += TraceSeed(a);
			if(N>=4)return 1;
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
		hits_in_seed.push_back(current_hit);
		if(hits_in_seed.size() >= MAX_SEED_HITS)break;
		Dtrk_hit *a = FindClosestXY(current_hit);
		if(!a)break;
		if(a->DistXY2(current_hit) > MAX_SEED_DIST2)break;
		if(last_hit){
			if(current_hit->CosPhiXY(a,last_hit) > 0.0)break;
		}
				
		last_hit = current_hit;
		current_hit = a;

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
	DQuickFit *fit = new DQuickFit();
	fit->SetMagneticFieldMap(bfield);
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *a = trkhits[i];
		if(!(a->flags&Dtrk_hit::IN_SEED))continue;
		fit->AddHitXYZ(a->x, a->y, a->z);
	}
	fit->FitCircle();
	x0 = fit->x0;
	y0 = fit->y0;
	r0 = sqrt(x0*x0 + y0*y0);
		
	// Set ON_CIRCLE flag for all hits close to the circle defined by x0,y0
	int N_in_seed_and_on_circle = 0;
	hits_on_circle.clear();
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *a = trkhits[i];
		a->flags &= ~(Dtrk_hit::ON_CIRCLE | Dtrk_hit::ON_LINE | Dtrk_hit::ON_TRACK);
		if(a->flags&Dtrk_hit::USED)continue;

		float dx = a->x - x0;
		float dy = a->y - y0;
		float r = sqrt(dx*dx + dy*dy);
			
		if(fabs(r0-r) <= MAX_CIRCLE_DIST){
			a->flags |= Dtrk_hit::ON_CIRCLE;
			hits_on_circle.push_back(a);
			if(a->flags & Dtrk_hit::IN_SEED)N_in_seed_and_on_circle++;
		}
	}
	
	// Record diagnostic (patfind) information
	if(dbg_seed_fit.size()<MAX_DEBUG_BUFFERS){
		dbg_in_seed.push_back(hits_in_seed);
		dbg_hoc.push_back(hits_on_circle);
		dbg_seed_fit.push_back(fit);
	}else{
		delete fit;
	}

	// If the IN_SEED hits don't actually end up on the circle
	// (I don't quite understand it, but it happens) then we should
	// set the ignore flag for the first IN_SEED hit and try again.
	if(N_in_seed_and_on_circle<3 || hits_on_circle.size()<4){
		ChopSeed();
		return 0;
	}

	return 1;
}

//------------------
// FindLineHits
//------------------
int DTrackCandidate_factory::FindLineHits(void)
{
	// The hits_on_circle vector should now contain pointers to
	// only those hits on the circle and there should be enough
	// hits (>=4) to do a fit
	
	// Find the phi-z angle and the z-vertex
	int ok = FindPhiZAngle();
	if(ok)ok = FindZvertex();	
	if(!ok){
		ChopSeed();
		return 0;
	}

	// Find all on-circle hits consistent with phizangle and zvertex
	float cos_phizangle = cos(phizangle);
	float sin_phizangle = sin(phizangle);
	hits_on_line.clear();
	//cout<<__FILE__<<":"<<__LINE__<<" phizangle="<<phizangle<<"  z_vertex="<<z_vertex<<" r0="<<r0<<endl;
	for(unsigned int i=0; i<hits_on_circle.size(); i++){
		Dtrk_hit *a = hits_on_circle[i];
		float dphi = a->phi_circle;
		float dz = a->z-z_vertex;
		float dr = sqrt(dphi*dphi + dz*dz);
		float sin_rel = (dphi*cos_phizangle - dz*sin_phizangle)/dr;
		float d = sin_rel*dr;
		//cout<<__FILE__<<":"<<__LINE__<<" d="<<d<<" z="<<a->z<<" dphi="<<dphi<<" dz="<<dz<<" dr="<<dr<<" sin_rel="<<sin_rel<<endl;
		if(fabs(d)<MAX_PHI_Z_DIST){
			// Flags hits as "ON_LINE" and push onto hits_on_line vector
			a->flags |= Dtrk_hit::ON_LINE;
			hits_on_line.push_back(a);
		}
	}
//cout<<__FILE__<<":"<<__LINE__<<" phizangle="<<phizangle<<" z_vertex="<<z_vertex<<endl;	
//cout<<__FILE__<<":"<<__LINE__<<" Nhits_on_line="<<hits_on_line.size()<<endl;

	if(hits_on_line.size()<4){
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
	Fill_phi_circle(hits_on_circle, x0, y0);

	// Loop over all hits on circle and fill the phizangle histo.
	phizangle_hist->Reset();
	for(unsigned int i=0; i<hits_on_circle.size(); i++){
		Dtrk_hit *a = hits_on_circle[i];

		float theta1 = atan2f((a->phi_circle - 0.0), a->z - TARGET_Z_MIN);
		float theta2 = atan2f((a->phi_circle - 0.0), a->z - TARGET_Z_MAX);
		float theta_min, theta_max;
		if(theta1<theta2){
			theta_min = theta1;
			theta_max = theta2;
		}else{
			theta_min = theta2;
			theta_max = theta1;
		}
		
		for(float t=theta_min; t<=theta_max; t+=phizangle_bin_size){
			phizangle_hist->Fill(t);
		}
	}
	
	// Use maximum bin for phi-z angle
	int xbin, ybin, zbin;
	phizangle_hist->GetMaximumBin(xbin,ybin,zbin);
	phizangle = phizangle_hist->GetXaxis()->GetBinCenter(xbin);
	
	// Bin content should be at least as musch as the minimum
	// number of hits required for a track.
	float height = phizangle_hist->GetBinContent(xbin);
	if(height<4)return 0;

	// For large angle tracks, the "peak" in the phizangle histo
	// becomes a plateau. In this case, the maximum bin returned
	// is the left-most bin of the plateau. (At least this  is the
	// current behaviour of ROOT.) Look for a plateau and use the
	// center of it if it exists.
	int Nbins = phizangle_hist->GetXaxis()->GetNbins();
	int xbin_right = xbin;
	for(int i=xbin+1; i<Nbins; i++){
		if(phizangle_hist->GetBinContent(i) != height)break;
		xbin_right = i;
	}
	phizangle_min = phizangle - phizangle_bin_size/2.0;
	phizangle_max = phizangle + phizangle_bin_size*(float)(xbin_right-xbin+1);
	phizangle = (phizangle_min + phizangle_max)/2.0;

	// Diagnostics
	if(dbg_phiz_hist.size()<MAX_DEBUG_BUFFERS){
		dbg_phiz_hist.push_back(new TH1F(*phizangle_hist));
		dbg_phiz_hist_seed.push_back(dbg_seed_fit.size()-1);
	}

	return 1;
}

//------------------
// Fill_phi_circle
//------------------
void DTrackCandidate_factory::Fill_phi_circle(vector<Dtrk_hit*> hits, float X, float Y)
{
	/// Fill in the phi_circle attribute for all Dtrk_hits in "hits" by
	/// calculating phi measured about the axis at x0,y0. Hits with either
	/// the IGNORE or USED flags set will be ignored.
	float r0 = sqrt(X*X + Y*Y);
	float x_last = -X;
	float y_last = -Y;
	float phi_last = 0.0;
	//unsigned int mask = Dtrk_hit::IGNORE | Dtrk_hit::USED;
//cout<<__FILE__<<":"<<__LINE__<<"  ------ x0="<<X<<" y0="<<Y<<" -----"<<endl;
	for(unsigned int i=0; i<hits.size(); i++){
		Dtrk_hit *a = hits[i];
		//if(a->flags & mask)continue;

		float dx = a->x - X;
		float dy = a->y - Y;
		float dphi = atan2f(dx*y_last - dy*x_last, dx*x_last + dy*y_last);
		float phi = phi_last +dphi;
		a->phi_circle = phi*r0;
//cout<<__FILE__<<":"<<__LINE__<<" z="<<a->z<<" dphi="<<a->phi_circle<<endl;
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
	zvertex_hist->Reset();

	// tan(x) has poles at x=+/- pi/2 What we really want
	// is cot(x) but that isn't always available. To avoid
	// the pole condition, calculate it using cos(x)/sin(x).
	float cot_phizangle_min = cos(phizangle_min)/sin(phizangle_min);
	float cot_phizangle_max = cos(phizangle_max)/sin(phizangle_max);

	for(unsigned int i=0; i<hits_on_circle.size(); i++){
		Dtrk_hit *a = hits_on_circle[i];

		// Find intersections with Z-axis for lines with min and
		// max phizangle going through point (a->z, a->phi_circle).
		float z1 = a->z - a->phi_circle*cot_phizangle_min;
		float z2 = a->z - a->phi_circle*cot_phizangle_max;
		float z_min, z_max;
		if(z1<z2){
			z_min = z1;
			z_max = z2;
		}else{
			z_min = z2;
			z_max = z1;
		}
		if(z_min<TARGET_Z_MIN) z_min = TARGET_Z_MIN;
		if(z_max>TARGET_Z_MAX) z_max = TARGET_Z_MAX;
		
		for(float z=z_min; z<=z_max; z+=z_vertex_bin_size){
			zvertex_hist->Fill(z);
		}
	}
	
	// Use maximum bin for z vertex
	int xbin, ybin, zbin;
	zvertex_hist->GetMaximumBin(xbin,ybin,zbin);
	z_vertex = zvertex_hist->GetXaxis()->GetBinCenter(xbin);
	
	// Bin content should be at least as musch as the minimum
	// number of hits required for a track.
	float height = zvertex_hist->GetBinContent(xbin);
	if(height<4)return 0;

	// Similar to phizangle method above, find range of plateau
	int Nbins = zvertex_hist->GetXaxis()->GetNbins();
	int xbin_right = xbin;
	for(int i=xbin+1; i<Nbins; i++){
		if(zvertex_hist->GetBinContent(i) != height)break;
		xbin_right = i;
	}
	float z_vertex_min = z_vertex - z_vertex_bin_size/2.0;
	float z_vertex_max = z_vertex + z_vertex_bin_size*(float)(xbin_right-xbin+1);
	z_vertex = (z_vertex_min + z_vertex_max)/2.0;

	// Diagnostics
	if(dbg_zvertex_hist.size()<MAX_DEBUG_BUFFERS){
		dbg_zvertex_hist.push_back(new TH1F(*zvertex_hist));
		dbg_zvertex_hist_seed.push_back(dbg_seed_fit.size()-1);
	}

	return 1;
}

//------------------
// FitTrack
//------------------
int DTrackCandidate_factory::FitTrack(void)
{
	// Create a new DTrackCandidate object and fill in the ihits
	// values as we loop over hits below.
	DTrackCandidate *trackcandidate = new DTrackCandidate;
	DQuickFit *fit = new DQuickFit();
	fit->SetMagneticFieldMap(bfield);
	for(unsigned int i=0; i<hits_on_line.size(); i++){
		Dtrk_hit *a = hits_on_line[i];
		
		// Add hit to fit
		fit->AddHitXYZ(a->x, a->y, a->z);
	}
	
	// Do a fit to the final subset of hits
	//fit->FitTrack();
	fit->FitTrack_FixedZvertex(z_vertex);
	
	trackcandidate->x0 = x0 = fit->x0;
	trackcandidate->y0 = y0 = fit->y0;
	r0 = sqrt(x0*x0 + y0*y0);
		
	// The following is from a fit of ratio of thrown to reconstructed
	// transverse momentum vs. theta
	double par[] = {0.984463, 0.150759, -0.414933, 0.257472, -0.055801};
	double theta = fit->theta;
	double ff = par[0]+theta*(par[1]+theta*(par[2]+theta*(par[3]+theta*par[4])));
	double sin_theta = sin(fit->theta);

	trackcandidate->z_vertex = fit->z_vertex;
	trackcandidate->dzdphi = r0/theta;
	trackcandidate->p_trans = fit->p_trans*ff;
	trackcandidate->p = trackcandidate->p_trans/sin_theta;
	trackcandidate->q = fit->q;
	trackcandidate->phi = fit->phi;
	trackcandidate->theta = fit->theta;

	// Do one last pass over the hits, marking all 
	// that are consistent with these parameters as used.
	MarkTrackHits(trackcandidate, fit);
	
	_data.push_back(trackcandidate);

	// For diagnostics
	if(dbg_track_fit.size()<MAX_DEBUG_BUFFERS){
		dbg_track_fit.push_back(fit);
		dbg_hol.push_back(hits_on_line);
		dbg_phizangle.push_back(phizangle);
		dbg_z_vertex.push_back(z_vertex);
		dbg_seed_index.push_back(dbg_seed_fit.size()-1);
	}else{
		delete fit;
	}

	return 1;
}

//------------------
// MarkTrackHits
//------------------
int DTrackCandidate_factory::MarkTrackHits(DTrackCandidate *trackcandidate, DQuickFit *fit)
{
	// The values of x0 and y0 are set in FitTrack as the
	// trackcandidate object is being filled. 
	Fill_phi_circle(trkhits, x0, y0);

	// Find all hits consistent with phizangle and zvertex
	float my_phizangle = -fit->q/fabs(fit->q)*tan(fit->theta);
//cout<<__FILE__<<":"<<__LINE__<<" my_phizangle/phizangle="<<my_phizangle/phizangle<<endl;
my_phizangle = phizangle;
	float cos_phizangle = cos(my_phizangle);
	float sin_phizangle = sin(my_phizangle);
	float r0_2pi_cos_phizangle = 2.0*M_PI*r0*cos_phizangle;
	hits_on_track.clear();
	int N_in_seed_and_on_track = 0;
	for(unsigned int i=0; i<trkhits.size(); i++){
		Dtrk_hit *a = trkhits[i];

		float x = a->x - x0;
		float y = a->y - y0;
		float r = sqrt(x*x + y*y);
		float dr = r - r0;
		
		// The phi angle could be off by an integral number
		// of 2pi's. To calculate the "real" distance, first
		// figure out how many 2pi's to shift by		
		float dphi = a->phi_circle;
		float dz = a->z-z_vertex;
		float d = dphi*cos_phizangle - dz*sin_phizangle;
		float n = floor(0.5 + d/r0_2pi_cos_phizangle);
		d -= n*r0_2pi_cos_phizangle;
//cout<<__FILE__<<":"<<__LINE__<<" z="<<a->z<<" d="<<d<<" d0="<<d+n*r0_2pi_cos_phizangle<<" dphi="<<dphi<<" n="<<d/r0_2pi_cos_phizangle<<" r0="<<r0<<endl;
		if(fabs(d)>MAX_PHI_Z_DIST)continue;
		if(fabs(dr)>MAX_CIRCLE_DIST)continue;
		
		// Flags hits as "ON_TRACK" and "USED" and push onto hits_on_track vector
		a->flags |= Dtrk_hit::ON_TRACK | Dtrk_hit::USED;
		hits_on_track.push_back(a);
		if(a->flags & Dtrk_hit::IN_SEED)N_in_seed_and_on_track++;

		// Add hit index to track in factory data
		trackcandidate->hitid.push_back(a->id);
	}

	// For diagnostics
	if(dbg_hot.size()<MAX_DEBUG_BUFFERS){
		dbg_hot.push_back(hits_on_track);
	}
	
	// It's possible that none of the seed hits are included in the
	// final track, when that happens, we have to ignore the first
	// seed hit on the next loop or else we'll keep getting the same
	// seed and the same result.
	if(hits_on_track.size()<3 || N_in_seed_and_on_track<1){
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

	printheader("        id: Nhits: x0(cm): y0(cm): z_vertex: dz/dphi:  q:     p: p_trans:   phi: theta:");

	for(unsigned int i=0; i<_data.size(); i++){
		DTrackCandidate *trackcandidate = _data[i];
		printnewrow();
		
		printcol("%d",    trackcandidate->id);
		printcol("%d",    trackcandidate->hitid.size());
		printcol("%3.1f", trackcandidate->x0);
		printcol("%3.1f", trackcandidate->y0);
		printcol("%3.1f", trackcandidate->z_vertex);
		printcol("%1.3f", trackcandidate->dzdphi);
		printcol("%+1.0f", trackcandidate->q);
		printcol("%3.3f", trackcandidate->p);
		printcol("%3.2f", trackcandidate->p_trans);
		printcol("%1.3f", trackcandidate->phi);
		printcol("%1.3f", trackcandidate->theta);

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
		cout<<i<<" x y z: "<<t->x<<" "<<t->y<<" "<<t->z<<endl;
		cout<<i<<" r phi: "<<t->r<<" "<<t->phi<<endl;
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


