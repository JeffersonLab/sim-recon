// $Id$
//
//    File: DFactory_DMCTrackCandidate_B.cc
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>

#include <TH1.h>

#include "DMCCheatHit.h"
#include "DFactory_DMCTrackCandidate_B.h"
#include "DQuickFit.h"



class TrkHitSort{
	public:
		bool operator()(DTrkHit* const &hit1, DTrkHit* const &hit2) const {
			return hit1->r > hit2->r;
		}
};
class TrkHitZSort{
	public:
		bool operator()(DTrkHit* const &hit1, DTrkHit* const &hit2) const {
			return hit1->z < hit2->z;
		}
};

//------------------
// DFactory_DMCTrackCandidate_B
//------------------
DFactory_DMCTrackCandidate_B::DFactory_DMCTrackCandidate_B()
{
	MAX_SEED_DIST = 5.0;
	MAX_SEED_DIST2 = MAX_SEED_DIST*MAX_SEED_DIST;
	MAX_CIRCLE_DIST = 2.0;
	MAX_PHI_Z_DIST = 10.0;
	MAX_DEBUG_BUFFERS = 0.0;
	TARGET_Z_MIN = 64.0;
	TARGET_Z_MAX = 66.0;
	
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
// fini
//------------------
derror_t DFactory_DMCTrackCandidate_B::fini(void)
{
	delete phizangle_hist;
	delete phi_relative;
	delete zvertex_hist;
	
	return NOERROR;
}

//------------------
// evnt
//------------------
derror_t DFactory_DMCTrackCandidate_B::evnt(DEventLoop *loop, int eventnumber)
{
	// Clear previous event from internal buffers
	ClearEvent();

	// Get the hits into the trkhits vector
	GetTrkHits(loop);

	// Loop as long as we keep finding tracks
	for(int i=0;i<10;i++){
		// Find a seed (group of hits that appear to be a clean track segment)
		if(!FindSeed())break;
		
		// Fit the seed hits to a circle and flag all unused hits
		// which fall close to the circle. If the fit fails, the
		// first IN_SEED hit automatically has it's IGNORE flag set
		// so we can just jump to the next iteration of the loop.
		if(!FitSeed())continue;

		// Using results of X/Y fit, find phi-z angle and z_vertex
		// Make a list of hits consistent with the values found.
		if(!FindLineHits())continue;
		
		// Fit the hits in the list created by above. Use result
		// to make a new DMCTrackCandidate
		if(!FitTrack())continue;
	}

	return NOERROR;
}

//------------------
// ClearEvent
//------------------
void DFactory_DMCTrackCandidate_B::ClearEvent(void)
{
	// Delete TrkHits from previous event
	for(unsigned int i=0; i<trkhits.size(); i++)delete trkhits[i];
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
void DFactory_DMCTrackCandidate_B::GetTrkHits(DEventLoop *loop)
{
	// Copy Cheat Hits into trkhits vector
	vector<const DMCCheatHit*> mccheathits;
	loop->Get(mccheathits);
	for(unsigned int i=0; i<mccheathits.size(); i++){
		const DMCCheatHit *mchit = mccheathits[i];
		if(!mchit->primary)continue;
		if(mchit->system>2)continue;
		float x = mchit->r*cos(mchit->phi);
		float y = mchit->r*sin(mchit->phi);
		float z = mchit->z;
		float r = mchit->r;
		float phi = mchit->phi;
		DTrkHit *hit = new DTrkHit(x,y,z,r,phi,mchit->system,mchit->track, i);
		trkhits.push_back(hit);
	}
	
	// Sort hits by r in X/Y plane
	trkhits_r_sorted = trkhits;
	sort(trkhits_r_sorted.begin(), trkhits_r_sorted.end(), TrkHitSort());

	// Order the track hits by z.
	sort(trkhits.begin(), trkhits.end(), TrkHitZSort());

}

//------------------
// FindSeed
//------------------
int DFactory_DMCTrackCandidate_B::FindSeed(void)
{

	// Loop over all un-used, non-ignored hits and try
	// tracing a seed from each. Once a seed with 4 or
	// more hits is found, return a 1 so a track can
	// be searched for.
	for(unsigned int i=0; i<trkhits_r_sorted.size(); i++){
		DTrkHit *a = trkhits_r_sorted[i];
		if(!(a->flags & (DTrkHit::USED | DTrkHit::IGNORE))){
		
			// Clear IN_SEED bit flag all hits
			for(unsigned int j=0; j<trkhits.size(); j++){
				trkhits[j]->flags &= ~(DTrkHit::IN_SEED);
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
int DFactory_DMCTrackCandidate_B::TraceSeed(DTrkHit *hit)
{
	// Starting with "hit", look for nearest neighbors to
	// trace the seed out on one side as far as possible
	// until we find a hit that is either too far away, or
	// changes directions by too large an angle
	int N = 0;
	DTrkHit *current_hit = hit;
	DTrkHit *last_hit = NULL;
	do{
		N++;
		current_hit->flags |= DTrkHit::IN_SEED;
		hits_in_seed.push_back(current_hit);
		DTrkHit *a = FindClosestXY(current_hit);
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
DTrkHit* DFactory_DMCTrackCandidate_B::FindClosestXY(DTrkHit *hit)
{
	DTrkHit *closest_hit = NULL;
	unsigned int mask = DTrkHit::USED | DTrkHit::IN_SEED | DTrkHit::IGNORE;
	float d2_min = 1000.0;
	for(unsigned int i=0; i<trkhits.size(); i++){
		DTrkHit *a = trkhits[i];
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
int DFactory_DMCTrackCandidate_B::FitSeed(void)
{
	// Do a quick circle fit to find the center of the seed
	DQuickFit *fit = new DQuickFit();
	for(unsigned int i=0; i<trkhits.size(); i++){
		DTrkHit *a = trkhits[i];
		if(!(a->flags&DTrkHit::IN_SEED))continue;
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
		DTrkHit *a = trkhits[i];
		a->flags &= ~(DTrkHit::ON_CIRCLE | DTrkHit::ON_LINE | DTrkHit::ON_TRACK);
		if(a->flags&DTrkHit::USED)continue;

		float dx = a->x - x0;
		float dy = a->y - y0;
		float r = sqrt(dx*dx + dy*dy);
			
		if(fabs(r0-r) <= MAX_CIRCLE_DIST){
			a->flags |= DTrkHit::ON_CIRCLE;
			hits_on_circle.push_back(a);
			if(a->flags & DTrkHit::IN_SEED)N_in_seed_and_on_circle++;
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
//cout<<__FILE__<<":"<<__LINE__<<" N_in_seed_and_on_circle="<<N_in_seed_and_on_circle<<"  hits_on_circle.size()"<<hits_on_circle.size()<<endl;
	if(N_in_seed_and_on_circle<3 || hits_on_circle.size()<4){
		ChopSeed();
		return 0;
	}

	return 1;
}

//------------------
// FindLineHits
//------------------
int DFactory_DMCTrackCandidate_B::FindLineHits(void)
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
		DTrkHit *a = hits_on_circle[i];
		float dphi = a->phi_circle;
		float dz = a->z-z_vertex;
		float dr = sqrt(dphi*dphi + dz*dz);
		float sin_rel = (dphi*cos_phizangle - dz*sin_phizangle)/dr;
		float d = sin_rel*dr;
		//cout<<__FILE__<<":"<<__LINE__<<" d="<<d<<" z="<<a->z<<" dphi="<<dphi<<" dz="<<dz<<" dr="<<dr<<" sin_rel="<<sin_rel<<endl;
		if(fabs(d)<MAX_PHI_Z_DIST){
			// Flags hits as "ON_LINE" and push onto hits_on_line vector
			a->flags |= DTrkHit::ON_LINE;
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
int DFactory_DMCTrackCandidate_B::FindPhiZAngle(void)
{
	// Fill the phi_circle field for all hits on the circle
	// centered at x0,y0
	Fill_phi_circle(hits_on_circle, x0, y0);

	// Loop over all hits on circle and fill the phizangle histo.
	phizangle_hist->Reset();
	for(unsigned int i=0; i<hits_on_circle.size(); i++){
		DTrkHit *a = hits_on_circle[i];

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
void DFactory_DMCTrackCandidate_B::Fill_phi_circle(vector<DTrkHit*> hits, float X, float Y)
{
	/// Fill in the phi_circle attribute for all DTrkHits in "hits" by
	/// calculating phi measured about the axis at x0,y0. Hits with either
	/// the IGNORE or USED flags set will be ignored.
	float r0 = sqrt(X*X + Y*Y);
	float x_last = -X;
	float y_last = -Y;
	float phi_last = 0.0;
	unsigned int mask = DTrkHit::IGNORE | DTrkHit::USED;
	for(unsigned int i=0; i<hits.size(); i++){
		DTrkHit *a = hits[i];
		if(a->flags & mask)continue;

		float dx = a->x - X;
		float dy = a->y - Y;
		float dphi = atan2f(dx*y_last - dy*x_last, dx*x_last + dy*y_last);
		float phi = phi_last +dphi;
		a->phi_circle = phi*r0;
		x_last = dx;
		y_last = dy;
		phi_last = phi;
	}
}

//------------------
// FindZvertex
//------------------
int DFactory_DMCTrackCandidate_B::FindZvertex(void)
{
	// Reset z_vertex histogram
	zvertex_hist->Reset();

	// tan(x) has poles at x=+/- pi/2 What we really want
	// is cot(x) but that isn't always available. To avoid
	// the pole condition, calculate it using cos(x)/sin(x).
	float cot_phizangle_min = cos(phizangle_min)/sin(phizangle_min);
	float cot_phizangle_max = cos(phizangle_max)/sin(phizangle_max);

	for(unsigned int i=0; i<hits_on_circle.size(); i++){
		DTrkHit *a = hits_on_circle[i];

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
		if(phizangle_hist->GetBinContent(i) != height)break;
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
int DFactory_DMCTrackCandidate_B::FitTrack(void)
{
	// Create a new DMCTrackCandidate object and fill in the ihits
	// values as we loop over hits below.
	DMCTrackCandidate *mctrackcandidate = new DMCTrackCandidate;
	DQuickFit *fit = new DQuickFit();
	for(unsigned int i=0; i<hits_on_line.size(); i++){
		DTrkHit *a = hits_on_line[i];
		
		// Add hit to fit
		fit->AddHitXYZ(a->x, a->y, a->z);
	}
	
	// Do a fit to the final subset of hits
	fit->FitTrack();
	
	mctrackcandidate->x0 = x0 = fit->x0;
	mctrackcandidate->y0 = y0 = fit->y0;
	r0 = sqrt(x0*x0 + y0*y0);
		
	mctrackcandidate->z_vertex = fit->z_vertex;
	mctrackcandidate->dphidz = fit->theta/r0;
	mctrackcandidate->p = fit->p;
	mctrackcandidate->p_trans = fit->p_trans;
	mctrackcandidate->q = fit->q;
	mctrackcandidate->phi = fit->phi;
	mctrackcandidate->theta = fit->theta;

	// Do one last pass over the hits, marking all 
	// that are consistent with these parameters as used.
	MarkTrackHits(mctrackcandidate);
	
	_data.push_back(mctrackcandidate);

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
int DFactory_DMCTrackCandidate_B::MarkTrackHits(DMCTrackCandidate *mctrackcandidate)
{
	// The values of x0 and y0 are set in FitTrack as the
	// mctrackcandidate object is being filled. 
	Fill_phi_circle(trkhits, x0, y0);

	// Find all hits consistent with phizangle and zvertex
	float cos_phizangle = cos(phizangle);
	float sin_phizangle = sin(phizangle);
	float r0 = sqrt(x0*x0 + y0*y0);
	hits_on_track.clear();
	int track_hits[10]={0,0,0,0,0,0,0,0,0,0};
	mctrackcandidate->Nhits = 0;
	for(unsigned int i=0; i<trkhits.size(); i++){
		DTrkHit *a = trkhits[i];

		float x = a->x - x0;
		float y = a->y - y0;
		float r = sqrt(x*x + y*y);
		float dr = r - r0;
		if(fabs(dr)>MAX_CIRCLE_DIST)continue;
		
		// The phi angle could be off by an integral number
		// of 2pi's. To calculate the "real" distance, first
		// figure out how many 2pi's to shift by		
		float dphi = a->phi_circle;
		float dz = a->z-z_vertex;
		float d = dphi*cos_phizangle - dz*sin_phizangle;
		float n = floor(0.5 + d/2.0/M_PI/r0);
		d -= n*2.0*M_PI*r0;
		if(fabs(d)>MAX_PHI_Z_DIST)continue;
		
		// Flags hits as "ON_TRACK" and "USED" and push onto hits_on_track vector
		a->flags |= DTrkHit::ON_TRACK | DTrkHit::USED;
		hits_on_track.push_back(a);
		if(a->track>=1 && a->track<=9)track_hits[a->track]++;

		// Add hit index to track in factory data
		mctrackcandidate->ihit[mctrackcandidate->Nhits++] = a->ihit;
		if(mctrackcandidate->Nhits>=MAX_IHITS){
			cout<<__FILE__<<":"<<__LINE__<<" More than "<<MAX_IHITS<<" hits on track. Truncating..."<<endl;
			break;
		}
	}
	
	// Fill in track and Ntrackhits fields of DMCTrackCandidate
	int imax = 0;
	for(int i=1;i<10;i++){
		if(track_hits[i]>track_hits[imax])imax = i;
	}
	mctrackcandidate->track = imax;
	mctrackcandidate->Ntrackhits = track_hits[imax];

	// For diagnostics
	if(dbg_hot.size()<MAX_DEBUG_BUFFERS){
		dbg_hot.push_back(hits_on_track);
	}
	
	if(hits_on_track.size()<3){
		ChopSeed();
		return 0;
	}

	return 1;
}

//------------------
// toString
//------------------
const string DFactory_DMCTrackCandidate_B::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: Nhits: x0(cm): y0(cm): z_vertex: dphi/dz:  q:   p: p_trans:   phi: theta:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DMCTrackCandidate *trackcandidate = _data[i];

		printnewrow();
		
		printcol("%d",    i);
		printcol("%d",    trackcandidate->Nhits);
		printcol("%3.1f", trackcandidate->x0);
		printcol("%3.1f", trackcandidate->y0);
		printcol("%3.1f", trackcandidate->z_vertex);
		printcol("%1.3f", trackcandidate->dphidz);
		printcol("%+1.0f", trackcandidate->q);
		printcol("%3.1f", trackcandidate->p);
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
// DebugMessage
//------------------
void DFactory_DMCTrackCandidate_B::DebugMessage(int line)
{
	// Find how many seed hits are from the same track
	int Nhits = 0;
	int in_seed=0, used=0, on_circle=0, ignore=0;
	int oc_hist[5]={0,0,0,0,0};
	int is_hist[5]={0,0,0,0,0};
	for(unsigned int i=0; i<trkhits.size(); i++){
		const DTrkHit *a = trkhits[i];
		Nhits++;
		if(a->flags&DTrkHit::IN_SEED)in_seed++;
		if(a->flags&DTrkHit::USED)used++;
		if(a->flags&DTrkHit::ON_CIRCLE)on_circle++;
		if(a->flags&DTrkHit::IGNORE)ignore++;
		
		if(a->track<1 || a->track>4)continue;
		if(a->flags&DTrkHit::IN_SEED)is_hist[a->track]++;
		if(a->flags&DTrkHit::ON_CIRCLE)oc_hist[a->track]++;
	}
	
	int oc_hits=0, oc_max_content=0;
	int is_hits=0, is_max_content=0;
	for(int i=1; i<=4; i++){
		oc_hits += oc_hist[i];
		is_hits += is_hist[i];
		if(oc_hist[i] > oc_max_content)oc_max_content = oc_hist[i];
		if(is_hist[i] > is_max_content)is_max_content = is_hist[i];
	}

	cout<<__FILE__<<":"<<line;
	cout<<" IN_SEED(good,bad)="<<in_seed<<"("<<is_max_content<<","<<is_hits-is_max_content<<")";
	cout<<" ON_CIRCLE(good,bad)="<<on_circle<<"("<<oc_max_content<<","<<oc_hits-oc_max_content<<")";
	cout<<" USED="<<used<<" IGNORE="<<ignore;
	cout<<" unused="<<trkhits.size()-used;
	cout<<endl;
}

//------------------
// SeedTrack
//------------------
int DFactory_DMCTrackCandidate_B::SeedTrack(void)
{
	// Find the track which most of the seed hits come from
	int is_hist[5]={0,0,0,0,0};
	for(unsigned int i=0; i<trkhits.size(); i++){
		const DTrkHit *a = trkhits[i];		
		if(a->track<1 || a->track>4)continue;
		if(a->flags&DTrkHit::IN_SEED)is_hist[a->track]++;
	}
	
	int track = 1;
	for(int i=2; i<=4; i++){
		if(is_hist[i] > is_hist[track])track = i;
	}
	return track;
}


