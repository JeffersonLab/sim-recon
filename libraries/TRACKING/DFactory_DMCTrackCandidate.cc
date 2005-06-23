// $Id$
//
//    File: DFactory_DMCTrackCandidate.cc
// Created: Sun Apr  3 12:38:16 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <string>
using namespace std;

#include <pthread.h>
#include <TThread.h>

#include "DFactory_DMCTrackCandidate.h"
#include "DMCCheatHit.h"
#include "DEventLoop.h"
#include "DQuickFit.h"
#include "DArcHit.h"

//------------------------------------------------------------------
// DFactory_DMCTrackCandidates (constructor)
//------------------------------------------------------------------
DFactory_DMCTrackCandidate::DFactory_DMCTrackCandidate()
{
	// Initialize storage vectors
	archits.clear();
	circles.clear();
	markers.clear();
	qfits.clear();
	
	// ROOT objects need unique names so they don't conflict with 
	// those defined by other factories in other threads.
	char str[256];
	unsigned int id = (unsigned int)pthread_self();

	// set limits for plot. This represents the space where the center 
	// of the circle can be. It can be (and often is) outside of the
	// bounds of the solenoid.
	circle_max = 150.0; // in cm.
	
	// The number of cm per bin (in one dimension) for the density histogram
	cm_per_bin = 4.5;

	// Creating and destroying histograms (and presumably other ROOT
	// objects) is not a thread safe operation (see chapter 21 pg 4
	// of the ROOT user's manual). Luckily, ROOT provides a locking
	// mechanism to prevent collisions.
	TThread::Lock();

	int Nbins = (int)floor(0.5 + 2.0*circle_max/cm_per_bin);
	sprintf(str,"density-%x", id);
	TH2F *density = new TH2F(str,"Density",Nbins,-circle_max,circle_max,Nbins,-circle_max,circle_max);	
	density_histos.push_back(density);

	// Histos of x and y coordinates of intersections of lines
	float intersect_circle_max = 500.0;
	Nbins = (int)floor(0.5 + 2.0*intersect_circle_max/cm_per_bin);
	sprintf(str,"density-x-%x", id);
	TH1F *density_x = new TH1F(str,"Intersection Density-X",Nbins,-intersect_circle_max,intersect_circle_max);	
	intersect_density_histos_x.push_back(density_x);
	sprintf(str,"density-y-%x", id);
	TH1F *density_y = new TH1F(str,"Intersection Density-Y",Nbins,-intersect_circle_max,intersect_circle_max);	
	intersect_density_histos_y.push_back(density_y);

	// max distance a line-of-circle-centers line can
	// be from a focal point and still be considered on the circle
	masksize = 5.0; // in cm
	masksize2 = masksize*masksize;
		
	// Create slope and intercept density histos
	sprintf(str,"slope-%x", id);
	TH1F *slope_density = new TH1F(str,"slope", 1000,-M_PI, +M_PI);
	slope_density_histos.push_back(slope_density);
	sprintf(str,"intercept-%x", id);
	TH1F *offset_density = new TH1F(str,"z intercept", 300, -100.0, +200.0);
	offset_density_histos.push_back(offset_density);
	
	// Release ROOT mutex
	TThread::UnLock();
}

//------------------------------------------------------------------
// ~DFactory_DMCTrackCandidates (destructor)
//------------------------------------------------------------------
DFactory_DMCTrackCandidate::~DFactory_DMCTrackCandidate()
{
	ClearEvent();
	
	// The ClearEvent() method doesn't delete the zeroth element histograms
	TThread::Lock();
	delete density_histos[0];
	delete intersect_density_histos_x[0];
	delete intersect_density_histos_y[0];
	delete slope_density_histos[0];
	delete offset_density_histos[0];
	TThread::UnLock();

	density_histos.clear();
	slope_density_histos.clear();
	offset_density_histos.clear();
}


//------------------------------------------------------------------
// ClearEvent
//------------------------------------------------------------------
void DFactory_DMCTrackCandidate::ClearEvent(void)
{
	for(unsigned int i=0; i<archits.size(); i++)delete archits[i];
	archits.clear();
	for(unsigned int i=0; i<circles.size(); i++)delete circles[i];
	circles.clear();
	for(unsigned int i=0; i<markers.size(); i++)delete markers[i];
	markers.clear();

	TThread::Lock();
	for(unsigned int i=1; i<density_histos.size(); i++)delete density_histos[i];
	density_histos.erase(density_histos.begin()+1, density_histos.end());
	for(unsigned int i=1; i<intersect_density_histos_x.size(); i++)delete intersect_density_histos_x[i];
	intersect_density_histos_x.erase(intersect_density_histos_x.begin()+1, intersect_density_histos_x.end());
	for(unsigned int i=1; i<intersect_density_histos_y.size(); i++)delete intersect_density_histos_y[i];
	intersect_density_histos_y.erase(intersect_density_histos_y.begin()+1, intersect_density_histos_y.end());
	for(unsigned int i=1; i<slope_density_histos.size(); i++)delete slope_density_histos[i];
	slope_density_histos.erase(slope_density_histos.begin()+1, slope_density_histos.end());
	for(unsigned int i=1; i<offset_density_histos.size(); i++)delete offset_density_histos[i];
	offset_density_histos.erase(offset_density_histos.begin()+1, offset_density_histos.end());
	TThread::UnLock();

	for(unsigned int i=0; i<qfits.size(); i++)delete qfits[i];
	qfits.clear();
	for(unsigned int i=0; i<intersect_points.size(); i++)delete intersect_points[i];
	intersect_points.clear();
	
}

//------------------
// evnt
//------------------
derror_t DFactory_DMCTrackCandidate::evnt(DEventLoop *loop, int eventnumber)
{
	// Empty vectors used to store objects for the event
	ClearEvent();

	// Get MCCheatHits and loop over them copying them into
	// the archit objects array
	vector<const DMCCheatHit*> mccheathits;
	loop->Get(mccheathits);
	for(unsigned int i=0; i<mccheathits.size(); i++){
		const DMCCheatHit *mccheathit = mccheathits[i];
	
		if(mccheathit->system!=1 && mccheathit->system!=2)continue;
		float x = mccheathit->r*cos(mccheathit->phi);
		float y = mccheathit->r*sin(mccheathit->phi);

		// Create new DArcHit object to hold info from this hit
		DArcHit *archit = new DArcHit();
		archit->track = mccheathit->track;
		archit->ihit = i;
		archit->SetXYZ(x,y,mccheathit->z);
		archits.push_back(archit);
	}
	if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<" Narchits:"<<archits.size()<<endl;
	if(archits.size()<3)return NOERROR;
	
	// sort the archits by z (mccheathits are sorted by track number , then z)
	sort(archits.begin(), archits.end(), ArcSort<DArcHit*>()); // ArcSort definied in DArcHit.h
	
	// Find circle patterns. When a circle is found, FindTracks is
	// automatically called to find all tracks associated with the circle.
	FindCircles();
	if(debug_level>0)cout<<__FILE__<<":"<<__LINE__<<" Ncircles:"<<circles.size()<<endl;
	
	// At this point, we seem to often have tracks which were "found"
	// multiple times (i.e. they have essentially the same fitted
	// parameters). Someday, I'll have to track that down, but until
	// then, I'll just go through and weed out duplicates.
	if(_data.size()>1){
		for(unsigned int i=0;i<_data.size()-1; i++){
			DMCTrackCandidate *a = _data[i];
			int filtered = 0;
			for(unsigned int j=i+1;j<_data.size(); j++){
				DMCTrackCandidate *b = _data[j];
				if(a->q == b->q)
					if(fabs(a->p_trans - b->p_trans)<0.05)
						if(fabs(a->phi - b->phi)<0.009)
							if(fabs(a->theta - b->theta)<0.009){
								ThereCanBeOnlyOne(i,j);
								filtered = 1;
								break;
							}
			}
			if(filtered){
				i--;
			}
		}
	}
	if(debug_level>0)cout<<__FILE__<<":"<<__LINE__<<" Ntracks:"<<_data.size()<<" (after filter)"<<endl;
	
	return NOERROR;
}

//------------------------------------------------------------------
// FindCircles
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::FindCircles(void)
{
	/// Call either FindCirclesMaskSub or FindCirclesHitSub
	///
	/// There are two methods for finding the peaks in the density
	/// histogram. The FindCirclesMaskSub routine will simply
	/// zero all bins within masksize of the bin with the maximum
	/// value. This method is probably quickest since you don't have
	/// to refill the density histogram over and over. The drawback
	/// is that there can be areas where high density occurs
	/// just from hits which have lines that are very close to the same
	/// slope so overlap in areas far from the focus. The result
	/// is false maxima being identified.
	///
	/// The second method is implemented in FindCirclesHitSub. This
	/// method regenerates the density histogram after finding
	/// each peak, but excludes hits assigned to tracks found from
	/// previous peaks. This is currently the method used and the
	/// only way to alter it is to modify the code.
	///
	/// Both of these methods may need to be optimized by adjusting the
	/// density function itself (DArcHit::Density()). At this point
	/// it seems as though it will be easier to look at contrasting
	/// the two if they are called from here.
	
	return FindCirclesIntersections();
	
	//return FindCirclesHitSub();
}

//------------------------------------------------------------------
// FindCirclesIntersections
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::FindCirclesIntersections(void)
{
	/// Find circle centers using 1-D histograms of coordinates of
	/// intersection points.

	// Clear the "used" flags on all archits
	for(unsigned int i=0;i<archits.size();i++ )archits[i]->used = 0;
	
	// Fill the intersect_points vector
	FindIntersectionPoints();

	// Copy pointer of first (default) density histos to local variables
	TH1F *density_x = intersect_density_histos_x[0];
	TH1F *density_y = intersect_density_histos_y[0];

	// Loop until we run out of circles
	float half_band_width = 3.0*masksize;
	int Nbins = density_y->GetXaxis()->GetNbins();
	int refill_primary_histo = 1;
	int band_direction = -1;
	TH1F *density_primary = NULL, *density_secondary = NULL;
	do{
		// If we did not find any tracks in the previous iteration then
		// we shouldn't re-fill the density histos since no more tracks
		// were marked as used. The previous peak is zeroed out at the
		// end of the previous iteration in th0se cases case.
		int xbin, ybin;
		if(refill_primary_histo){
			// Fill histos with values from intersection points of unused hits
			density_x->Reset();
			density_y->Reset();
			for(unsigned int i=0; i<intersect_points.size(); i++){
				intersect_pt_t *pt = intersect_points[i]; // is this efficient?
				if(pt->archit_a->used)continue;
				if(pt->archit_b->used)continue;
				density_x->Fill(pt->x);
				density_y->Fill(pt->y);
			}
			
			// This is just for debugging
			if((int)intersect_density_histos_x.size()<=max_density_histos){
				intersect_density_histos_x.push_back(new TH1F(*density_x));
			}
			if((int)intersect_density_histos_y.size()<=max_density_histos){
				intersect_density_histos_y.push_back(new TH1F(*density_y));
			}
			
			// Here we choose whether we need to refill the x density with
			// a cut on y or the other way around. We choose by whichever
			// projection's maximum bin is smallest. That projection is
			// more likely to have the resonances spread out so we're less likely
			// to have multiple tracks overlapping.
			int tmp1,tmp2;
			density_x->GetMaximumBin(xbin,tmp1,tmp2);
			density_y->GetMaximumBin(ybin,tmp1,tmp2);
			if(density_x->GetBinContent(xbin) < density_y->GetBinContent(ybin)){
				band_direction = BAND_DIR_X;
				density_primary = density_x;
				density_secondary = density_y;
				if(debug_level>10)cout<<__FILE__<<":"<<__LINE__<<" Band direction: X"<<endl;
			}else{
				band_direction = BAND_DIR_Y;
				density_primary = density_y;
				density_secondary = density_x;
				if(debug_level>10)cout<<__FILE__<<":"<<__LINE__<<" Band direction: Y"<<endl;
			}
		}
		
		// Center of peak in band direction 
		int prim_bin,tmp1,tmp2;
		density_primary->GetMaximumBin(prim_bin,tmp1,tmp2);
		if(density_primary->GetBinContent(prim_bin) < 10)break;
		float prim_center = density_primary->GetXaxis()->GetBinCenter(prim_bin);
		
		// Fill the secondary histo with values from un-used hits which have
		// values within +/- 3 masksizes of peak in the primary direction
		density_secondary->Reset();
		for(unsigned int i=0; i<intersect_points.size(); i++){
			intersect_pt_t *pt = intersect_points[i]; // is this efficient?
			if(pt->archit_a->used)continue;
			if(pt->archit_b->used)continue;
			float vprim = band_direction==BAND_DIR_X ? pt->x:pt->y;
			if(vprim > prim_center+half_band_width)continue;
			if(vprim < prim_center-half_band_width)continue;
			density_secondary->Fill(band_direction==BAND_DIR_X ? pt->y:pt->x);
		}
		

		// Find all peaks in the secondary density histo, zeroing out peaks
		// as we fail to find tracks in them. If a track is found, we jump
		// back to the outer loop and refill the x/y histos and choose a new
		// primary.
		int Ntracks = 0;
		do{
			// Center of peak in direction orthogonal to band
			int sec_bin,tmp1,tmp2;
			density_secondary->GetMaximumBin(sec_bin,tmp1,tmp2);
			if(density_secondary->GetBinContent(sec_bin) < 10)break;
			float sec_center = density_secondary->GetXaxis()->GetBinCenter(sec_bin);

			// Assign X/Y depending upon what the band direction is
			float x,y;
			if(band_direction==BAND_DIR_X){
				x = prim_center;
				y = sec_center;
			}else{
				y = prim_center;
				x = sec_center;
			}
			
			// Try this X/Y combination
			Ntracks = FindTrack_RoughXY(x,y);
			
			// If a track was found, break here and let the density_x
			// histogram be re-filled on the next iteration of the outer
			// loop, excluding the points used by the recently found track(s)
			if(Ntracks > 0)break;
			
			// Zero out the density_secondary histo for 5 bins on either side
			// of the peak so we can try again with the next peak.
			for(int i=sec_bin-5;i<=sec_bin+5;i++){
				if(i<1 || i>Nbins)continue;
				density_secondary->SetBinContent(i, 0.0);
			}
			
			// At some point, the entire histogram will get zeroed out
			// and we will break the loop at the top
		}while(1);
		
		// If no tracks were found then no more hits were marked as new
		// and we'll infinitely loop unless something is changed. If no
		// tracks were found, we zero out the current peak in the primary
		// and clear the refill flag so it is not refilled at the next
		// iteration.
		if(Ntracks==0){
			refill_primary_histo = 0;
			for(int i=prim_bin-5;i<=prim_bin+5;i++){
				if(i<1 || i>Nbins)continue;
				density_primary->SetBinContent(i, 0.0);
			}
		}else{
			refill_primary_histo = 1;
		}
		
		// Similar to the secondary-loop above, the entire histogram will 
		// eventually get zeroed out and we will break the loop at the top
	}while(1);

	return NOERROR;
}

//------------------------------------------------------------------
// FindIntersectionPoints
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::FindIntersectionPoints(void)
{
	
	// This is redundant since it is also done in ClearEvent(). Oh well...
	intersect_points.clear();
	
	// Loop over all pairs of lines, finding and recording the intersection
	// point of each.
	for(unsigned int i=0;i<archits.size()-1;i++){
		DArcHit *a = archits[i];
		for(unsigned int j=i+1;j<archits.size();j++){
			DArcHit *b = archits[j];
	
			// Intersection of two lines:
			// c1*x + c2*y = c3
			// d1*x + d2*y = d3
			intersect_pt_t *pt = new intersect_pt_t;
			float c1 = a->orientation==DArcHit::Y_OF_X ? -a->m:1.0;
			float c2 = a->orientation==DArcHit::Y_OF_X ? 1.0:-a->m;
			float c3 = a->b;
			float d1 = b->orientation==DArcHit::Y_OF_X ? -b->m:1.0;
			float d2 = b->orientation==DArcHit::Y_OF_X ? 1.0:-b->m;
			float d3 = b->b;
			pt->x = (d2*c3 - d3*c2)/(d2*c1 - c2*d1);
			pt->y = (d3*c1 - d1*c3)/(d2*c1 - c2*d1);
			
			if(finite(pt->x) && finite(pt->y)){
				pt->archit_a = a;
				pt->archit_b = b;
				intersect_points.push_back(pt);
			}
		}
	}

	return NOERROR;
}

//------------------------------------------------------------------
// FindTracks_RoughXY
//------------------------------------------------------------------
int DFactory_DMCTrackCandidate::FindTrack_RoughXY(float x, float y)
{
	/// This routine takes as input a rough guess of the X/Y coordinates of
	/// the center of a circle. It fits a circle to the points whose
	/// corresponding archit line comes close to the "rough" center
	/// in order to find an accurate center. It then sets the on_circle
	/// flag along with the delta_phi attribute for all un-used archits
	/// with lines that come close to the
	/// accurate center. Finally, it calls FindTracks().


	// The maxmimum bin is not terribly accurate as the center of the
	// circle. Use DQuickFit to quickly find a better center.
	DQuickFit *fit = new DQuickFit();
	for(unsigned int i=0; i<archits.size(); i++){
		DArcHit *a = archits[i];
		if(a->used)continue;
		float d2 = a->Dist2ToLine(x,y);
		if(d2 <= masksize2){
			fit->AddHit(a->rhit, a->phihit, a->zhit);
		}
	}
		
	if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<" NHits for circle fit:"<<fit->GetNhits()<<endl;
	int Ntracks = 0;
	if(fit->GetNhits()>=3){
		fit->FitCircle();
		float x0 = -fit->x0;	// why do we need the minus sign?
		float y0 = -fit->y0;	// why do we need the minus sign?
		float phi0 = atan2(y0, x0);
		if(phi0<0.0)phi0 += 2.0*M_PI;

		// Record the location of the maximum
		TEllipse *circle = new TEllipse();
		circle->SetX1(x0);
		circle->SetY1(y0);
		circle->SetR1(masksize);
		circle->SetR2(masksize);
		circles.push_back(circle);
		
		// Loop over all hits. Use the "on_circle" flag in the DArchit objects to
		// mark hits which actually appear to be on this circle.
		int Nhits_on_circle = 0;
		float delta_phi_offset = 0.0;
		float last_delta_phi=-1000.0;
		float r0 = sqrt(x0*x0 + y0*y0);
		for(unsigned int i=0;i<archits.size();i++){
			DArcHit *a = archits[i];
			if(a->used){
				a->on_circle = 0;
				continue;
			}
			float d2 = a->Dist2ToLine(x0,y0);
			if(d2 > masksize2){
				a->on_circle = 0;
				continue;
			}
			a->on_circle = 1;
			Nhits_on_circle++;

			// Find relative angle between this hit
			// and vector pointing to origin
			float delta_x = a->xhit+x0;
			float delta_y = a->yhit+y0;
		
			// Calculate delta_phi and force it to be in the 0 to +2PI range
			float delta_phi = atan2(delta_y, delta_x) - phi0;
			while(delta_phi<0.0)delta_phi += 2.0*M_PI;
			if(delta_phi>M_PI)delta_phi -= 2.0*M_PI;
			if(last_delta_phi != -1000.0){
				float dphi = delta_phi - last_delta_phi;
				if(fabs(dphi) > M_PI){
					delta_phi_offset += dphi<0.0 ? +2.0*M_PI:-2.0*M_PI;
					if(debug_level>20)cout<<__FILE__<<":"<<__LINE__<<"  delta_phi_offset="<<delta_phi_offset<<endl;
				}
			}
			last_delta_phi = delta_phi;
			delta_phi += delta_phi_offset;
			a->delta_phi = delta_phi*r0; // put in cm so it has same units as a->zhit
		}
		if(debug_level>10)cout<<__FILE__<<":"<<__LINE__<<" Nhits on circle #"<<circles.size()<<": "<<Nhits_on_circle<<endl;

		// Now look for tracks in phi/z plane using this circle as the
		// axis for phi.
		Ntracks = FindTrack(x0,y0);
	}
	delete fit;

	return Ntracks;
}

//------------------------------------------------------------------
// FindTracks
//------------------------------------------------------------------
int DFactory_DMCTrackCandidate::FindTrack(float x0, float y0)
{
	// At this point the archits which fall on the circle centered at
	// x0,y0 should all have their "on_circle" flags set. Now we need
	// to identify those hits which fall on a line in the phi/z plane
	// and fit them to a 3D track.
	// 
	// There are several approaches we can take here. All of them
	// involve histograming the slope and intercept of all possible
	// pairs of phi/z points from them on_circle hits. One major
	// choice to be made here is whether or not to try looping over
	// peaks in one or both of the histograms while looking for a
	// track. Note that we are already looping over peaks in the
	// intersection density histos in FindCirclesIntercetion(). Ideally,
	// we would loop over all peaks in all density histos while looking
	// for tracks, but in the order in which we are most likely to find
	// them. In other words, test the secondary peaks in
	// the slope histo only after checking the secondary peaks in the
	// x and y intersection histos.
	//
	// Looping over peaks in the slope and density histos here gave
	// poorer results at reduced efficiency from just simulataneously
	// filling the slope and offset histos and taking the largest
	// peak in each. Hence, we stick with that.
	if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"   Finding tracks in circle at x="<<x0<<" y="<<y0<<endl;

	// Copy pointer of first (default) density histos to local variables
	TH1F *slope_density = slope_density_histos[0];
	TH1F *offset_density = offset_density_histos[0];
	
	// Fill the offset histo
	FillOffsetDensityHisto();

	// Find the largest peak in the offset density 
	int xbin, ybin, zbin;
	offset_density->GetMaximumBin(xbin,ybin,zbin);
	if(debug_level>10)cout<<__FILE__<<":"<<__LINE__<<"   offset density: maxbin="<<offset_density->GetXaxis()->GetBinCenter(xbin)<<"  bin contents="<<offset_density->GetBinContent(xbin)<<endl;
	if(offset_density->GetBinContent(xbin) < 3)return 0;
	float z_vertex = offset_density->GetXaxis()->GetBinCenter(xbin);
z_vertex=65.0;
	if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"    z_vertex = "<<z_vertex<<endl;
	if(z_vertex<-95.0 || z_vertex>+195.0)return 0;
	
	// Fill slope density histo using the found z_vertex
	FillSlopeDensityHisto(z_vertex);

	// Find the largest peak in the slope density 
	slope_density->GetMaximumBin(xbin,ybin,zbin);
	if(debug_level>10)cout<<__FILE__<<":"<<__LINE__<<"   slope density: maxbin="<<slope_density->GetXaxis()->GetBinCenter(xbin)<<"  bin contents="<<slope_density->GetBinContent(xbin)<<endl;
	if(slope_density->GetBinContent(xbin) < 3)return 0;
	float phi_z_angle = slope_density->GetXaxis()->GetBinCenter(xbin);
	if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"    slope_density = "<<phi_z_angle<<endl;

	// Now try and make the track
	int Ntracks = MakeTrack(phi_z_angle, z_vertex, sqrt(x0*x0+y0*y0));
	
	// If we got here, no tracks were found
	return Ntracks;
}

//------------------------------------------------------------------
// FillOffsetDensityHisto
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::FillOffsetDensityHisto(void)
{
	/// Calculate parameters for all possible pairs of phi-z coordinates
	/// for all archits currently flagged as being "on_circle" recording
	/// the results in the phi_z_lines vector.
	
	// Copy pointer of first (default) density histos to local variables
	TH1F *offset_density = offset_density_histos[0];
	
	offset_density->Reset();
	for(unsigned int i=0; i<archits.size()-1; i++){
		DArcHit *a = archits[i];
		if(!a->on_circle)continue;
		for(unsigned int k=i+1; k<archits.size(); k++){ 
			DArcHit *b = archits[k];
			if(!b->on_circle)continue;
			
			// archits are sorted by Z so "a" is always upstream of "b"
			float delta_z = b->zhit - a->zhit;
			float delta_phi = b->delta_phi-a->delta_phi; //NOTE: already scaled by r0 (see FindTrack_RoughXY)
			float phi_z_angle = atan2(delta_phi, delta_z);
			float z_vertex = a->zhit - a->delta_phi*delta_z/delta_phi;

			//if(debug_level>10)cout<<__FILE__<<":"<<__LINE__<<"  phi_z_angle="<<phi_z_angle<<"  z_vertex="<<z_vertex<<endl;
			if(!finite(phi_z_angle) || !finite(z_vertex))continue;
			
			// Limit target range (Probably shouldn't hardcode this
			if(fabs(z_vertex - 65.0)>30.0)continue;
			
			offset_density->Fill(z_vertex);
		}
	}
	
	// This is just for debugging
	if((int)offset_density_histos.size()<=max_density_histos){
		offset_density_histos.push_back(new TH1F(*offset_density));
	}

	return NOERROR;
}

//------------------------------------------------------------------
// FillSlopeDensityHisto
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::FillSlopeDensityHisto(float z_vertex)
{
	/// Calculate the slope between each on_circle archit and the
	/// specified z_vertex position on the beamline filling the slope
	/// density histogram.
	
	// Copy pointer of first (default) density histos to local variables
	TH1F *slope_density = slope_density_histos[0];
	
	slope_density->Reset();
	for(unsigned int i=0; i<archits.size()-1; i++){
		DArcHit *a = archits[i];
		if(!a->on_circle)continue;
			
		// angle between this hit and the vertex
		float delta_z = a->zhit - z_vertex;
		float delta_phi = a->delta_phi; //NOTE: already scaled by r0 (see FindTrack_RoughXY)
		float phi_z_angle = atan2(delta_phi, delta_z);

		if(!finite(phi_z_angle))continue;
						
		slope_density->Fill(phi_z_angle);
	}
	
	// This is just for debugging
	if((int)slope_density_histos.size()<=max_density_histos){
		slope_density_histos.push_back(new TH1F(*slope_density));
	}

	return NOERROR;
}
//------------------------------------------------------------------
// MakeTrack
//------------------------------------------------------------------
int DFactory_DMCTrackCandidate::MakeTrack(float phi_z_angle, float z_vertex, float r0)
{
		
	// Create a new DQuickFit object to do the fit (should there
	// be enough hits to do one).
	DQuickFit *fit = new DQuickFit();
		
	// Create a new DMCTrackCandidate object and fill in the ihits
	// values as we loop over hits below. If it
	// turns out we don't have enough hits after looping over them
	// all, then we will delete the object.
	DMCTrackCandidate *mctrackcandidate = new DMCTrackCandidate;
	mctrackcandidate->Nhits = 0;

	// At this point, the "on_circle" flag in the array of DArcHits should
	// indicate which hits are consistent with the current circle
	// (i.e. centered at x0,y0). We now look for all hits with the "on_circle"
	// flag set that also have a delta_phi and zhit within a
	// certain distance of the line defined by phi_z_angle and z_vertex.

	for(unsigned int i=0;i<archits.size();i++){
		DArcHit *a = archits[i];
		if(!a->on_circle)continue;
		if(a->used)continue;

		// Calculate distance of the archit phi/z point to the line
		// defined by phi_z_angle and z_vertex.
		float cos_phi_z_angle = cos(phi_z_angle);
		float d = a->delta_phi*cos_phi_z_angle - (a->zhit-z_vertex)*sin(phi_z_angle);
		if(fabs(d)>masksize){
			// sometimes, the delta_phi values can be off by a multiple
			// of 2pi. Try and recover those points.
			float fN = round(d/(r0*2.0*M_PI*cos_phi_z_angle));
			d -= fN*2.0*M_PI*r0*cos_phi_z_angle;
		}
		if(debug_level>10)cout<<__FILE__<<":"<<__LINE__<<"        distance["<<i<<"] = "<<d<<"   N 2pis="<<d/(r0*2.0*M_PI*cos(phi_z_angle))<<endl;
		if(fabs(d)>masksize)continue; // this point is too far away from the line
		
		// Add hit to DQuickFit object
		fit->AddHit(a->rhit, a->phihit, a->zhit);
					
		// Add hit index to track in factory data
		mctrackcandidate->ihit[mctrackcandidate->Nhits++] = a->ihit;
				
		// Flag this hit as having been used
		a->used = 1;
				
		if(mctrackcandidate->Nhits>=MAX_IHITS){
			cout<<__FILE__<<":"<<__LINE__<<" More than "<<MAX_IHITS<<" hits on track. Truncating..."<<endl;
			break;
		}
	}			
		
	// If enough hits were added, then do the fit and record
	// the results
	int NTracks = 0;
	if(fit->GetNhits()>=3){
		fit->FitTrack();
		mctrackcandidate->x0 = fit->x0;
		mctrackcandidate->y0 = fit->y0;
		float r0 = sqrt(fit->x0*fit->x0 + fit->y0*fit->y0);
		
		mctrackcandidate->z_vertex = fit->z_vertex;
		mctrackcandidate->dphidz = fit->theta/r0;
		mctrackcandidate->p = fit->p;
		mctrackcandidate->p_trans = fit->p_trans;
		mctrackcandidate->q = fit->q;
		mctrackcandidate->phi = fit->phi;
		mctrackcandidate->theta = fit->theta;
		_data.push_back(mctrackcandidate);
		NTracks++;
		if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"      Adding track (Nhits="<<fit->GetNhits()<<")  -"<<_data.size()<<"-"<<endl;

		// Keep the DQuickFit object around
		qfits.push_back(fit);
	}else{
		// Oops! not enough hits for a track. Delete this one.
		delete mctrackcandidate;
		delete fit;
		if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"      Aborting track (Nhits="<<fit->GetNhits()<<")"<<endl;
	}

	return NTracks;
}

//------------------------------------------------------------------
// GetQFit
//------------------------------------------------------------------
DQuickFit* DFactory_DMCTrackCandidate::GetQFit(int n)
{
	if(n<0 || n>=(int)qfits.size())return NULL;
	
	return qfits[n];
}

//------------------------------------------------------------------
// DrawPhiZPoints
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::DrawPhiZPoints(int which)
{
	/// This is for debugging/development only.
	///
	/// Draw markers on the current canvas. The coordinates
	/// are the phi value on the Y-axis and the z value of the
	/// hit on the X-axis. The value of phi is relative to the
	/// focus. All hits passing within masksize of a focus are
	/// plotted in a color corresponding to the focus. A single hit
	/// can thus be plotted more than once, but will necessarily
	/// show up in different places on the plot for the different
	/// foci because phi will be different.
	///
	/// The idea here is that points who are really from the same
	/// track will fall on a line. This just helps visualize
	/// this for development. (The program patfind uses this).
	
	for(unsigned int i=0; i<markers.size(); i++)delete markers[i];
	markers.clear();

	int colors[] = {kRed,kBlue,kMagenta,kGreen,kBlack};
	for(unsigned int j=0;j<circles.size();j++){
		if(which>0 && (int)j!=which-1)continue;
		TEllipse *circle = circles[j];
		float x0 = circle->GetX1();
		float y0 = circle->GetY1();
		float phi0 = atan2(y0, x0);
		if(phi0<0.0)phi0 += 2.0*M_PI;
		float delta_phi_offset = 0.0;
		float last_delta_phi=-1000.0;
		for(unsigned int i=0;i<archits.size();i++){
			DArcHit *a = archits[i];
			
			if(a->Dist2ToLine(x0,y0) > masksize2)continue;
			
			// Find relative angle between this hit
			// and vector pointing to origin
			float delta_x = a->xhit+x0;
			float delta_y = a->yhit+y0;
		
			// Calculate delta_phi and force it to be in the 0 to +2PI range
			float delta_phi = atan2(delta_y, delta_x) - phi0;
			while(delta_phi<0.0)delta_phi += 2.0*M_PI;
			if(delta_phi>M_PI)delta_phi -= 2.0*M_PI;
			if(last_delta_phi != -1000.0){
				float dphi = delta_phi - last_delta_phi;
				if(fabs(dphi) > M_PI){
					delta_phi_offset += dphi<0.0 ? +2.0*M_PI:-2.0*M_PI;
					if(debug_level>20)cout<<__FILE__<<":"<<__LINE__<<"  delta_phi_offset="<<delta_phi_offset<<endl;
				}
			}
			last_delta_phi = delta_phi;
			delta_phi += delta_phi_offset;
			
			TMarker *marker = new TMarker();
			marker->SetX(a->zhit);
			marker->SetY(delta_phi);
			marker->SetMarkerColor(colors[j%6]);
			marker->SetMarkerStyle(20+j);
			marker->Draw();
			markers.push_back(marker);

			if(markers.size()>=300)break;
		}
		if(markers.size()>=300)break;
	}

	return NOERROR;
}

//------------------------------------------------------------------
// SetNumDensityHistograms
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::SetMaxDensityHistograms(int N)
{
	/// Set the number of density histograms to keep.
	///
	/// This is meant as a diagnostic and should not normally be used.
	/// It tells the factory to copy and keep around multiple 
	/// histograms corresponding to subsequent stages of the track
	/// finding. This can cost a little in memory, but the real
	/// cost is in the CPU time to copy the histograms. This value
	/// is actually used by the intersection density X , density Y
	/// slope, and z_vertex histos so changing it can have a large affect.

	max_density_histos = N;

	return NOERROR;
}

//------------------------------------------------------------------
// GetDensityHistogram
//------------------------------------------------------------------
TH2F* DFactory_DMCTrackCandidate::GetDensityHistogram(int n)
{
	/// Return a pointer to 2-D density histogram n where n is value from
	/// 0 to 1 less than the last call to SetNumDensityHistograms().
	
	if(n<0 || n>=(int)density_histos.size())return NULL;

	return density_histos[n];
}

//------------------------------------------------------------------
// GetIntersectDensityHistogramX
//------------------------------------------------------------------
TH1F* DFactory_DMCTrackCandidate::GetIntersectDensityHistogramX(int n)
{
	/// Return a pointer to 1-D density histogram representing X 
	/// distribution of intersection points
	
	if(n<0 || n>=(int)intersect_density_histos_x.size())return NULL;

	return intersect_density_histos_x[n];
}

//------------------------------------------------------------------
// GetIntersectDensityHistogramY
//------------------------------------------------------------------
TH1F* DFactory_DMCTrackCandidate::GetIntersectDensityHistogramY(int n)
{
	/// Return a pointer to 1-D density histogram representing Y
	/// distribution of intersection points
	
	if(n<0 || n>=(int)intersect_density_histos_y.size())return NULL;

	return intersect_density_histos_y[n];
}

//------------------------------------------------------------------
// GetSlopeDensityHistogram
//------------------------------------------------------------------
TH1F* DFactory_DMCTrackCandidate::GetSlopeDensityHistogram(int n)
{
	/// Return a pointer to 1-D slope density histogram n where n is value from
	/// 0 to 1 less than the last call to SetNumDensityHistograms().
	
	if(n<0 || n>=(int)slope_density_histos.size())return NULL;

	return slope_density_histos[n];
}

//------------------------------------------------------------------
// GetOffsetDensityHistogram
//------------------------------------------------------------------
TH1F* DFactory_DMCTrackCandidate::GetOffsetDensityHistogram(int n)
{
	/// Return a pointer to 1-D z-offset density histogram n where n is value from
	/// 0 to 1 less than the last call to SetNumDensityHistograms().
	
	if(n<0 || n>=(int)offset_density_histos.size())return NULL;

	return offset_density_histos[n];
}

//------------------------------------------------------------------
// ThereCanBeOnlyOne
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::ThereCanBeOnlyOne(int trk1, int trk2)
{
	/// See the comment at the end of evnt(). Basically, we need to choose
	/// which one track to keep and adjust the _data array to keep it.
	
	// If we're keeping the track further down the list, then copy
	// its results into the position near the front of the list.
	// We should be using the chisq from the fits, but that isn't
	// being calculated for the 3D track at the moment.
	DQuickFit *qf1 = qfits[trk1];
	DQuickFit *qf2 = qfits[trk2];
	DQuickFit *qf;
	int trk = -1;
	if(qf1->GetNhits() < qf2->GetNhits()){
		trk = trk2;
		qf = qf2;
	}else{
		trk = trk1;
		qf = qf1;
	}
		
	// Delete the losing track and shift the ones after it up in the list(s)
	delete _data[trk];
	_data.erase(_data.begin() + trk);
	delete qf;
	qfits.erase(qfits.begin() + trk);

	return NOERROR;	
}

//------------------
// toString
//------------------
const string DFactory_DMCTrackCandidate::toString(void)
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


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// ---------- Things below are unused -----------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
#if 0

//------------------------------------------------------------------
// FindCirclesHitSub
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::FindCirclesHitSub(void)
{
	/// Find circles by repeatedly filling the density
	/// histogram and looking for the maximum. As each
	/// maximum is found, record its position in the
	/// circles[] array and look for tracks in the phi/z
	/// plane (by calling FindTracks(x,y)). Hits used in
	/// 3-D tracks are flagged as used and do not contribute
	/// to filling the density histo on subsequent iterations.

	// Clear the "used" flags on all archits
	for(unsigned int i=0;i<archits.size();i++ )archits[i]->used = 0;

	// Loop until we run out of circles
	do{
		// Copy pointer to first (default) density histo to local variable
		TH2F *density_histo = density_histos[0];
	
		// Fill the density histogram
		FillArcDensityHistogram(density_histo);
		if((int)density_histos.size()<=max_density_histos)density_histos.push_back(new TH2F(*density_histo));

		// Find the coordinates of the maximum
		int xbin, ybin, zbin;
		density_histo->GetMaximumBin(xbin,ybin,zbin);
		float x = density_histo->GetXaxis()->GetBinCenter(xbin);
		float y = density_histo->GetYaxis()->GetBinCenter(ybin);
		
		// Use our "rough" x/y values to find a better set an then
		// look for tracks
		int Ntracks = FindTracks_RoughXY(x, y);
		
		// If no tracks were found, then no more hits were marked "used"
		// meaning we will loop infinitely (if the loop is not broken out
		// of below). Set the used flag for all hits within masksize
		// of this circle if no tracks were found. (There may be a something
		// more clever to do here, but right now I'm not sure what).
		if(Ntracks<1){
			for(unsigned int i=0;i<archits.size();i++){
				DArcHit *a = archits[i];
				if(a->used)continue;
				float d2 = a->Dist2ToLine(x,y);
				if(d2 <= masksize2)a->used = 1;
			}
		}
		
		// Count how many un-used hits we have left
		int Nhits_not_used = 0;
		for(unsigned int i=0;i<archits.size();i++){
			DArcHit *a = archits[i];
			if(!a->used)Nhits_not_used++;
		}
		if(debug_level>2)cout<<__FILE__<<":"<<__LINE__<<"   Number of unused hits: "<<Nhits_not_used<<endl;
					
		// If less than 4 hits remain unused, we can stop looking now
		if(Nhits_not_used<4)break;
		
		// I ran across a strange event in which 6 hits were unused, but
		// none of them came close to the circle center. This caused an
		// infinite loop. We break that possibility here.
		static int last_Nhits_not_used=0;
		if(Nhits_not_used == last_Nhits_not_used)break;
		last_Nhits_not_used = Nhits_not_used;

		// With 20 circles, something just ain't right and we should bail
	}while(circles.size()<20);
	
	return NOERROR;
}

//------------------------------------------------------------------
// FillArcDensityHistogram
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::FillArcDensityHistogram(TH2F *hist)
{
	/// Loop over all archits and fill the density histo for
	/// any that don't have their used flag set.
	hist->Reset();

	for(unsigned int i=0;i<archits.size();i++){
		DArcHit *a = archits[i];
		if(a->used)continue;
		
		a->FillArcDensityHistogram(hist);
		//break;
	}
	
	return NOERROR;
}

//------------------------------------------------------------------
// FindTracks
//------------------------------------------------------------------
int DFactory_DMCTrackCandidate::FindTracks(float x0, float y0)
{
	/// Find lines in the phi-z plane and create tracks from them.
	///
	/// This gets called from FindCircles (or FindCirclesHitSub)
	/// whenever a circle is found. The x,y coordinates given are
	/// used as the axis about which to calculate phi.
	/// A list of phi,z points is then
	/// made. The slope and z-intercept of every possible pair of
	/// phi,z points are then used to fill 2 1-dimensional histograms.
	/// Since a single circle in the x,y plane may correspond to
	/// multiple tracks, then multiple peaks may appear in the slope
	/// histogram (one for each track). If the tracks came from different
	/// vertices, then multiple peaks may also appear in the 
	/// z-intercept histogram. For now though, we assume only a 
	/// single vertex.
	///
	/// Each peak in the slope histo is used along with the peak in
	/// the z-intercept histogram to determine hits which form a track
	/// (i.e. ones close to the line defined by the slope/intercept.)
	/// When a group of hits are determined to form a track, they
	/// are used to create a DQuickFit object which is then used to
	/// quickly fit track parameters in 3-D. The results are used to
	/// fill the factory contents for this factory.

	if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"   Finding tracks in circle at x="<<x0<<" y="<<y0<<endl;

	// Fill in the phi_z_lines array
	FindPhiZLines();
	
	TH1F *slope_density = slope_density_histos[0];
	TH1F *offset_density = offset_density_histos[0];
	int Nslope_density_bins = slope_density->GetXaxis()->GetNbins();

	// Fill the slope and z-intercept histos for all hits consistent
	// with a circle at x0,y0.
	float r0 = sqrt(x0*x0 + y0*y0);
	FillSlopeIntDensityHistos(x0, y0);
		
	// For debugging
	if((int)slope_density_histos.size()<=max_density_histos){
		slope_density_histos.push_back(new TH1F(*slope_density));
		offset_density_histos.push_back(new TH1F(*offset_density));
	}

	// For now, assume a single vertex for all tracks included in
	// in this circle. Just use the maximum
	// in the z-intercept histo as the z coordinate of the vertex.
	int z_bin = offset_density_histos[0]->GetMaximumBin();
	float z_vertex = offset_density_histos[0]->GetBinCenter(z_bin);
	if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"    z_vertex from offset density histo: "<<z_vertex<<endl;

	// Find all peaks in the slope histogram
	int Ntracks_found = 0;
	do{
		if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"    Looking for track in SlopeInt histo"<<endl;
			
		// Use a simple algorithm here. Just use the maximum as the
		// peak center. Zero out the histogram near the peak to look
		// for a new peak. The number of pairs of phi,z points
		// is N(N-1)/2 . Thus, we need to find N which is the number
		// of unused archits on the circle.
		int Nuseable_archits = 0;		
		for(unsigned int i=0;i<archits.size();i++){
			DArcHit *a = archits[i];
			if(!a->on_circle || a->used)continue;
			Nuseable_archits++;
		}
		
		// Set the limit on the peak size. We do this so we can catch
		// the tracks with very few hits which all resonate.
		int min_peak_size = (Nuseable_archits*(Nuseable_archits-1))/2/2;
		if(min_peak_size > 20)min_peak_size = 20;
		if(min_peak_size < 3)min_peak_size = 3;
		
		int slope_bin = slope_density_histos[0]->GetMaximumBin();
		if(slope_density_histos[0]->GetBinContent(slope_bin)<min_peak_size){
			if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"      No more tracks slope_bin="<<slope_bin<<"  bin content="<<slope_density_histos[0]->GetBinContent(slope_bin)<<"  min_peak_size="<<min_peak_size<<endl;
			break;
		}
		if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"      Peak in slope histo: slope_bin="<<slope_bin<<"  bin content="<<slope_density_histos[0]->GetBinContent(slope_bin)<<"  min_peak_size="<<min_peak_size<<endl;
		
		// (see note in FillSlopeIntDensityHistos about 600/M_PI)
		float phi_z_angle = slope_density_histos[0]->GetBinCenter(slope_bin);
			
		// Create a new DQuickFit object to do the fit (should there
		// be enough hits to do one).
		DQuickFit *fit = new DQuickFit();
		
		// Create a new DMCTrackCandidate object and fill in the ihits
		// values as we loop over hits below. If it
		// turns out we don't have enough hits after looping over them
		// all, then we will delete the object.
		DMCTrackCandidate *mctrackcandidate = new DMCTrackCandidate;
		mctrackcandidate->Nhits = 0;

		// At this point, the "on_circle" flag in the array of DArcHits should
		// indicate which hits are consistent with the current circle
		// (i.e. centered at x0,y0). We now look for all hits with the "on_circle"
		// flag set that also have a delta_phi and zhit within a
		// certain distance of the line defined by slope and z_vertex.

		for(unsigned int i=0;i<archits.size();i++){
			DArcHit *a = archits[i];
			if(!a->on_circle)continue;

			// Calculate distance to line.
			// The value of phi_z_angle represents the angle of the
			// line in the phi-z plane make wrt to the beam line (x-axis).
			// Oddly enough, the distance of an arbitrary point to
			// this line is linear in the point's coordinates.
			// (see note in FillSlopeIntDensityHistos about 600/M_PI)
			float d = 600.0/M_PI*a->delta_phi*cos(phi_z_angle) - (a->zhit-z_vertex)*sin(phi_z_angle);

			//if(debug_level>10)cout<<__FILE__<<":"<<__LINE__<<"        masksize="<<masksize<<"    d="<<d<<"  ratio="<<((a->zhit-z_vertex)*sin(phi_z_angle))/(600.0/M_PI*a->delta_phi*cos(phi_z_angle))<<" sin_theta="<<sin(phi_z_angle)<<" cos_theta="<<cos(phi_z_angle)<<endl;

			if(fabs(d)<masksize){
				// Add hit to DQuickFit object
				fit->AddHit(a->rhit, a->phihit, a->zhit);
					
				// Add hit index to track in factory data
				mctrackcandidate->ihit[mctrackcandidate->Nhits++] = a->ihit;
				
				// Flag this hit as having been used
				a->used = 1;
				
				if(mctrackcandidate->Nhits>=MAX_IHITS){
					cout<<__FILE__<<":"<<__LINE__<<" More than "<<MAX_IHITS<<" hits on track. Truncating..."<<endl;
					break;
				}
			}
		}
			
		// If enough hits were added, then do the fit and record
		// the results
		if(fit->GetNhits()>=3){
			fit->FitTrack();
			mctrackcandidate->x0 = fit->x0;	// why do we need the minus sign?
			mctrackcandidate->y0 = fit->y0;	// why do we need the minus sign?
		
			mctrackcandidate->z_vertex = fit->z_vertex;
			mctrackcandidate->dphidz = fit->theta/r0;
			mctrackcandidate->p = fit->p;
			mctrackcandidate->p_trans = fit->p_trans;
			mctrackcandidate->q = fit->q;
			mctrackcandidate->phi = fit->phi;
			mctrackcandidate->theta = fit->theta;
			_data.push_back(mctrackcandidate);
			Ntracks_found++;
			if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"      Adding track (Nhits="<<fit->GetNhits()<<")  -"<<_data.size()<<"-"<<endl;

			// Keep the DQuickFit object around
			qfits.push_back(fit);
		}else{
			// Oops! not enough hits for a track. Delete this one.
			delete mctrackcandidate;
			delete fit;
			if(debug_level>1)cout<<__FILE__<<":"<<__LINE__<<"      Aborting track (Nhits="<<fit->GetNhits()<<")"<<endl;
		}
			

		// Zero out the slope_density histo for 5 bins on either side
		// of the peak.
		for(int i=slope_bin-5;i<=slope_bin+5;i++){
			if(i<1 || i>Nslope_density_bins)continue;
			slope_density->SetBinContent(i, 0.0);
		}
	}while(1);
	

	return Ntracks_found;
}

//------------------------------------------------------------------
// FillSlopeIntDensityHistos
//------------------------------------------------------------------
derror_t DFactory_DMCTrackCandidate::FillSlopeIntDensityHistos(float x0, float y0)
{
	/// This is called by FindTracks()
	
	TH1F *slope_density = slope_density_histos[0];
	TH1F *offset_density = offset_density_histos[0];
	slope_density->Reset();
	offset_density->Reset();
	float phi0 = atan2(y0, x0);
	if(phi0<0.0)phi0 += 2.0*M_PI;
	
	// Loop over all hits. Use the "on_circle" flag in the DArchit objects to
	// record who is and isn't used to fill the slope and z-intercept histos
	for(unsigned int i=0;i<archits.size();i++){
		DArcHit *a = archits[i];
		if(a->used){
			a->on_circle = 0;
			continue;
		}
		float d2 = a->Dist2ToLine(x0,y0);
		if(d2 > masksize2){
			a->on_circle = 0;
			continue;
		}
		a->on_circle = 1;

		// Find relative angle between this hit
		// and vector pointing to origin
		float x = a->xhit+x0;
		float y = a->yhit+y0;
		
		// Calculate delta_phi and force it to be in the -PI to +PI range
		a->delta_phi = atan2(y, x)-phi0;
		while(a->delta_phi<0.0)a->delta_phi += 2.0*M_PI;
		if(a->delta_phi>M_PI)a->delta_phi -= 2.0*M_PI;
	}

	// Loop over all pairs of phi,z points, filling the slope
	// and z-intercept histos
	
	for(unsigned int i=0; i<archits.size()-1; i++){
		DArcHit *a = archits[i];
		if(!a->on_circle)continue;
		for(unsigned int k=i+1; k<archits.size(); k++){ 
			DArcHit *b = archits[k];
			if(!b->on_circle)continue;
			if(b->delta_phi<0.0 && a->delta_phi>0.0)continue; // filter out pairs which can't fall on the same line
			if(b->delta_phi>0.0 && a->delta_phi<0.0)continue; // filter out pairs which can't fall on the same line

			// archits are sorted by Z so "a" is always upstream of "b"
			float delta_z = b->zhit - a->zhit;
			if(fabs(delta_z)==0.0)continue; // don't include hits from same plane
			float delta_phi = b->delta_phi-a->delta_phi;

			// OK, this may need to be re-thought, but here goes:
			// There is a problem using the slope due to the pole
			// at 90 degrees. Using the angle in phi vs. z space
			// is tricky since the typical range of values 
			// for phi is much smaller than that of z. This leads
			// to very poor sensitivity for forward going tracks.
			// We alleviate the problem by scaling
			// up the value of delta_phi to increase sensitivity 
			// in small phi-z angles. This is at the cost of becoming
			// less sensitive to large (~ 90 degree) angles.
			// This sensitivity comes about from using a histogram
			// of fixed bin size over the entire range of angles.
			
			float m = atan2(600.0/M_PI*(double)delta_phi, (double)delta_z);
			float z = a->zhit - a->delta_phi*delta_z/delta_phi;
			if(!finite(m) || !finite(z))continue;
			
			slope_density->Fill(m);
			offset_density->Fill(z);
		}
	}
	
	return NOERROR;
}

//------------------------------------------------------------------
// FindCirclesMaskSub
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FindCirclesMaskSub(void)
{
	/// Find circles by filling the density histogram and looking
	/// for the maximum. As each maximum is found, record its position in the
	/// circles[] array and zero out the masksize area around the
	/// maximum in the density histo.

	Ncircles = 0;	

	// Fill the density histogram
	TH2F *tmp = density_histos[0];
	TAxis *xaxis = tmp->GetXaxis();
	TAxis *yaxis = tmp->GetYaxis();
	int Nxbins = xaxis->GetNbins();
	int Nybins = yaxis->GetNbins();
	FillArcDensityHistogram(tmp);

	do{

		// Find the coordinates of the maximum
		int xbin, ybin, zbin;
		tmp->GetMaximumBin(xbin,ybin,zbin);
		float x = xaxis->GetBinCenter(xbin);
		float y = yaxis->GetBinCenter(ybin);
		
		// The maxmimum bin is not terribly accurate as the center of the
		// circle. Use DQuickFit to quickly find a better center.
		DQuickFit *fit = new DQuickFit();
		DArcHit *a = archit;
#if 0
		for(int i=0;i<Narchits;i++ ,a++){
			float d = a->DistToLine(x,y);
			if(d <= masksize){
				fit->AddHit(a->rhit, a->phihit, a->zhit);
			}
		}
		if(fit->GetNhits()>=3){
			fit->FitCircle();
			x = -fit->x0;	// why do we need the minus sign?
			y = -fit->y0;	// why do we need the minus sign?
		}
#endif
		delete fit;

		// Count the hits within masksize of the maximum
		int Nhits_this_track = 0;
		int Nhits_not_used = 0;
		a = archit;
		for(int i=0;i<Narchits;i++ ,a++){
			float d = a->DistToLine(x,y);
			if(d > masksize){
				if(!a->used)Nhits_not_used++;
				continue;
			}else{
				a->used = 1;
				Nhits_this_track++;
			}
		}
		
		// There should be at least 4 hits for it to be a track
		if(Nhits_this_track>=4){
		
			// Record the location of the maximum
			circles[Ncircles].SetX1(x);
			circles[Ncircles].SetY1(y);
			circles[Ncircles].SetR1(masksize);
			circles[Ncircles++].SetR2(masksize);
		}

		// Zero all bins within masksize of the maximum.
		float maxval = tmp->GetBinContent(xbin,ybin,zbin);
		cout<<__FILE__<<":"<<__LINE__<<" maxval = "<<maxval<<endl;
#if 0
		ZeroNeighbors(tmp, xbin,ybin);
#else
		int Nmaskbins = (int)ceil(masksize/cm_per_bin);
		for(int i=xbin-Nmaskbins; i<xbin+Nmaskbins; i++){
			if(i<1 || i>Nxbins)continue;
			for(int j=ybin-Nmaskbins; j<ybin+Nmaskbins; j++){
				if(j<1 || j>Nybins)continue;
				float d = sqrt(pow((double)(i-xbin),2.0) + pow((double)(j-ybin),2.0));
				if(d*cm_per_bin > masksize)continue;
				
				tmp->SetBinContent(i,j,0.0);
			}
		}
#endif
		
		if(maxval < 15.0)break;

	}while(Ncircles<20);
	
	return NOERROR;
}

//------------------------------------------------------------------
// ZeroNeighbors
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::ZeroNeighbors(TH2F *hist, int xbin, int ybin)
{
	/// Zero all of the 3x3 group of bins of in hist centered on
	/// xbin, ybin. Any bins above a set limit will call us again
	/// so this is a reentrant routine.
	int Nxbins = hist->GetXaxis()->GetNbins();
	int Nybins = hist->GetYaxis()->GetNbins();
	for(int i=xbin-1; i<xbin+1; i++){
		if(i<1 || i>Nxbins)continue;
		for(int j=ybin-2; j<ybin+2; j++){
			if(j<1 || j>Nybins)continue;
			
			float val = hist->GetBinContent(i,j,0);
			hist->SetBinContent(i,j,0.0);
			if(val > 10.0)ZeroNeighbors(hist,i,j);
		}
	}	

	return NOERROR;
}

//------------------------------------------------------------------
// FindCirclesInt
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FindCirclesInt(void)
{
	/// This essentially works the same as FindCirclesHitSub() only
	/// it avoids filling and dealing with the large 2-D histograms.
	/// It does this by looking at the intersection points of all
	/// possible pairs of lines and finding ones in which many other
	/// lines pass near.
	
	// Clear the "used" flag on all DArcHits
	DArcHit *a = archit;
	for(int i=0;i<Narchits;i++ ,a++)a->used=0;

	// Loop until all tracks are found
	Ncircles = 0;	
	do{
		
		// Loop over all possible intersection points (of the unused
		// hits) and find the one which has the most lines passing
		// within masksize.
		float x0, y0;
		int max_lines_within_masksize = 0;
		DArcHit *a = archit;
		for(int i=0;i<Narchits;i++ ,a++){
			if(a->used)continue;

			DArcHit *b = a;
			b++;
			for(int j=i+1;j<Narchits;j++ ,b++){
				if(b->used)continue;
				
				float x,y;
				int n = IntersectionDensity(a,b,x,y);
				if(n>max_lines_within_masksize){
					max_lines_within_masksize = n;
					x0 = x;
					y0 = y;
				}
			}
		}
		
		// Need a minimum number of hits to be considered a track
		if(max_lines_within_masksize<4)break;
		
		// OK, looks like we found a good intersection. Go through
		// the hits again and do a circle fit to find a better
		// x0, y0 value to use as the true focus.
		DQuickFit *fit = new DQuickFit();
		DArcHit *c = archit;
		for(int k=0;k<Narchits;k++ ,c++){
			if(c->used)continue;
			if(c->Dist2ToLine(x0,y0)>masksize2)continue;
			fit->AddHit(c->rhit, c->phihit, c->zhit);
		}
		fit->FitCircle();
		x0 = -fit->x0;	// why do we need the minus sign?
		y0 = -fit->y0;	// why do we need the minus sign?
		delete fit;

		// Loop over the hits one final time and mark all that
		// are within masksize of the focus as used
		c = archit;
		for(int k=0;k<Narchits;k++ ,c++){
			if(c->used)continue;
			if(c->Dist2ToLine(x0,y0)>masksize2)continue;
			c->used = 1;
		}
		
		// Finally, record the focus
		if(Ncircles<32){
			circles[Ncircles].SetX1(x0);
			circles[Ncircles].SetY1(y0);
			circles[Ncircles].SetR1(masksize);
			circles[Ncircles++].SetR2(masksize);
		}
	}while(Ncircles<32);	
	
	return NOERROR;
}

//------------------------------------------------------------------
// IntersectionDensity
//------------------------------------------------------------------
int DFactory_MCTrackCandidates::IntersectionDensity(DArcHit *a, DArcHit *b, float &x, float&y)
{
	// Intersection of two lines:
	// c1*x + c2*y = c3
	// d1*x + d2*y = d3
	float c1 = a->orientation==DArcHit::Y_OF_X ? -a->m:1.0;
	float c2 = a->orientation==DArcHit::Y_OF_X ? 1.0:-a->m;
	float c3 = a->b;
	float d1 = b->orientation==DArcHit::Y_OF_X ? -b->m:1.0;
	float d2 = b->orientation==DArcHit::Y_OF_X ? 1.0:-b->m;
	float d3 = b->b;
	x = (d2*c3 - d3*c2)/(d2*c1 - c2*d1);
	y = (d3*c1 - d1*c3)/(d2*c1 - c2*d1);
		
	// It's possible that x or y could be infinite (parallel lines)
	// In this case, just skip to the next hit.
	if(!finite(x) || !finite(y))return 0;

	DArcHit *c = archit;
	int lines_within_masksize = 2;
	for(int k=0;k<Narchits;k++ ,c++){
		if(c->used)continue;
		if(c==a || c==b)continue;
		if(c->Dist2ToLine(x,y)<=masksize2)lines_within_masksize++;
	}

	return lines_within_masksize;
}
#endif // 0
