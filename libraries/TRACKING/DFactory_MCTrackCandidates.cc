// $Id$

#include "DEvent.h"
#include "DFactory_MCTrackCandidates.h"
#include "DQuickFit.h"
#include "DFactory_MCCheatHits.h"

//------------------------------------------------------------------
// DFactory_MCTrackCandidates (constructor)
//------------------------------------------------------------------
DFactory_MCTrackCandidates::DFactory_MCTrackCandidates(DEvent *event):DFactory(event, "MCTrackCandidates", sizeof(MCTrackCandidate_t))
{
	// Initialize Nlines Nellipse;
	Narchits = 0;
	Ncircles = 0;
	Nmarkers = 0;
	markers = NULL;
	
	// set limits for plot. This represents the space where the center 
	// of the circle can be. It can be (and often is) outside of the
	// bounds of the solenoid.
	circle_max = 150.0; // in cm.
	
	// The number of bins per cm (in one dimension) for the density histogram
	bins_per_cm = 4.0;

	int Nbins = (int)round(2.0*circle_max/bins_per_cm);
	density = new TH2F("density","Density",Nbins,-circle_max,circle_max,Nbins,-circle_max,circle_max);	
	Ndensity_histos = 1;
	density_histos[0] = density;

	// max distance a line-of-circle-centers line can
	// be from a focal point and still be considered on the circle
	masksize = 2.0; // in cm
	
	// See note in FindCircles
	flip_x_axis = 0;
	
	// Create slope and intercept density histos
	slope_density = new TH1F("slope","slope", 1000,-0.05,0.05);
	slope_density_histos[0] = new TH1F(*slope_density);
	offset_density = new TH1F("intercept","z intercept", 2100, -100.0,2000.0);
	offset_density_histos[0] = new TH1F(*offset_density);
}

//------------------------------------------------------------------
// ~DFactory_MCTrackCandidates (destructor)
//------------------------------------------------------------------
DFactory_MCTrackCandidates::~DFactory_MCTrackCandidates()
{
	for(int i=0; i<Ndensity_histos; i++)delete density_histos[i];
	delete slope_density;
	delete offset_density;
	if(Nmarkers)delete markers;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::evnt(int eventnumber)
{	
	// Get MCCheatHits and loop over them copying them into
	// the archit objects array
	DContainer *mccheathits = event->Get("MCCheatHits");
	MCCheatHit_t *mccheathit = (MCCheatHit_t*)mccheathits->first();
	Narchits = 0;
	for(int i=0;i<mccheathits->nrows;i++, mccheathit++){
		if(Narchits>=300)break;
	
		if(mccheathit->system!=1 && mccheathit->system!=2)continue;
		float x = mccheathit->r*cos(mccheathit->phi);
		float y = mccheathit->r*sin(mccheathit->phi);
		
		// When drawing these circles on the screen, the natural coordinates
		// of the screen have the x-axis pointing to the right. However, in
		// the lab coordinate system, the x-axis points to the left (when
		// looking downstream). If the flip_x_axis flag is set, change the
		// sign of the x-coordinate so it is displayed properly.
		// (This was put in here originally for patfind, but may be useful
		// for other event viewers)
		if(flip_x_axis)x = -x;
		
		archit[Narchits].track = mccheathit->track;
		archit[Narchits].ihit = i;
		archit[Narchits++].SetXYZ(x,y,mccheathit->z);
	}
	
	// Find circle patterns first. (The results are left in circles[])
	FindCircles();
	
	// Split tracks found in X/Y using the Z-info
	FillSlopeIntDensityHistos();
	
	// Fill in Factory info
	// THIS IS NOT VERY EFFICIENT RIGHT NOW!! It re-calculates the
	// distance from the line-of-circle-centers to the focus which
	// is done in at least 2 other places before we get here. See 
	// FindCirclesHitSub(), FindCirclesMaskSub(), and FillSlopeIntDensityHistos()
	for(int j=0;j<Ncircles;j++){
		MCTrackCandidate_t *mctrackcandidate = (MCTrackCandidate_t*)_data->Add();
		
		float x0 = circles[j].GetX1();
		float y0 = circles[j].GetY1();

		mctrackcandidate->x0 = x0;
		mctrackcandidate->y0 = y0;
		
		// Not implemented yet
		mctrackcandidate->z_vertex = 0.0;
		mctrackcandidate->dphidz = 0.0;

		DArcHit *a = archit;
		mctrackcandidate->Nhits = 0;
		for(int i=0;i<Narchits;i++ ,a++){
			
			float d = a->DistToLine(x0,y0);
			if(d > masksize)continue;
			
			mctrackcandidate->ihit[mctrackcandidate->Nhits++] = a->ihit;
		}
	}	

	return NOERROR;
}

//------------------------------------------------------------------
// FindCircles
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FindCircles(void)
{
	/// Call either FindCirclesMaskSub or FindCirclesHitSub
	///
	/// There are two methods for finding the peaks in the density
	/// histogram. The FindCirclesMaskSub routine will simply
	/// zero all bins within masksize of the bin with the maximum
	/// value. This method is probably quickest since you don't have
	/// to refill the desnity histogram over and over. The drawback
	/// is that there can be areas where high density occurs
	/// just from hits which have lines that are very close to the same
	/// slope so overlap in areas far from the focus. The result
	/// is false maxima being identified.
	///
	/// The second method is implemented in FindCirclesHitSub. This
	/// method regenerates the density histogram after finding
	/// each peak, but excludes hits from contributing if their 
	/// lines passed within masksize of the focus. The drawback
	/// of this method is that for events in which one track
	/// has hits whose lines happen to pass over the focus of another
	/// hit, they will be removed and the peak at their focus will 
	/// disappear before being identified.
	///
	/// Both of these methods should be optimized by adjusting the
	/// density fundtion itself (DArcHit::Density()). At this point
	/// it seems as though it will be easier to look at contrasting
	/// the two if they are called from here.
	
	return FindCirclesHitSub();
}

//------------------------------------------------------------------
// FindCirclesHitSub
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FindCirclesHitSub(void)
{
	/// Find circles by repeatedly filling the density
	/// histogram and looking for the maximum. As each
	/// maximum is found, record its position in the
	/// circles[] array and remove the contributing hits
	/// from the density histo.

	Ncircles = 0;	
	do{
		// use tmp so we can more easily build in a mechanism for
		// filling multiple histograms for debugging
		TH2F *tmp = density_histos[Ncircles<Ndensity_histos ? Ncircles:Ndensity_histos-1];
	
		// Fill the density histogram
		FillArcDensityHistogram(tmp);

		// Find the coordinates of the maximum
		int xbin, ybin, zbin;
		tmp->GetMaximumBin(xbin,ybin,zbin);
		float x = tmp->GetXaxis()->GetBinCenter(xbin);
		float y = tmp->GetYaxis()->GetBinCenter(ybin);
		
		// The maxmimum bin is not terribly accurate as the center of the
		// circle. Use DQuickFit to quickly find a better center.
		DQuickFit *fit = new DQuickFit();
		DArcHit *a = archit;
		for(int i=0;i<Narchits;i++ ,a++){
			float d = a->DistToLine(x,y);
			if(d <= 2.0*masksize){
				fit->AddHit(a->rhit, a->phihit, a->zhit);
			}
		}
		if(fit->GetNhits()>=3){
			fit->FitCircle();
			x = -fit->x0;	// why do we need the minus sign?
			y = -fit->y0;	// why do we need the minus sign?
		}
		delete fit;

		// Flag all hits within masksize of the circle center as used
		int Nhits_this_track = 0;
		int Nhits_not_used = 0;
		a = archit;
		for(int i=0;i<Narchits;i++ ,a++){
			if(a->used)continue;
			float d = a->DistToLine(x,y);
			if(d <= masksize){
				a->used = 1;
				Nhits_this_track++;
			}else{
				Nhits_not_used++;
			}
		}
		
		// If there are less than 4 hits on this track, assume that
		// means it's not really a track and we've already found them
		// all.
		if(Nhits_this_track<4)break;
		
		// Record the location of the maximum
		circles[Ncircles].SetX1(x);
		circles[Ncircles].SetY1(y);
		circles[Ncircles].SetR1(masksize);
		circles[Ncircles++].SetR2(masksize);
					
		// If less than 4 hits remain unused, we can stop looking now
		if(Nhits_not_used<4)break;

	}while(Ncircles<20);
	
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
		int Nmaskbins = (int)ceil(masksize*bins_per_cm);
		for(int i=xbin-Nmaskbins; i<xbin+Nmaskbins; i++){
			if(i<1 || i>Nxbins)continue;
			for(int j=ybin-Nmaskbins; j<ybin+Nmaskbins; j++){
				if(j<1 || j>Nybins)continue;
				float d = sqrt(pow((double)(i-xbin),2.0) + pow((double)(j-ybin),2.0));
				if(d/bins_per_cm > masksize)continue;
				
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
		float masksize2 = masksize*masksize;
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
	float masksize2 = masksize*masksize;
	for(int k=0;k<Narchits;k++ ,c++){
		if(c->used)continue;
		if(c==a || c==b)continue;
		if(c->Dist2ToLine(x,y)<=masksize2)lines_within_masksize++;
	}

	return lines_within_masksize;
}

//------------------------------------------------------------------
// FillArcDensityHistogram
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FillArcDensityHistogram(TH2F *hist)
{
	/// Loop over all archits and fill the density histo for
	/// any that don't have their used flag set.
	hist->Reset();
	DArcHit *a = archit;
	for(int i=0;i<Narchits;i++ ,a++){
		if(a->used)continue;
		
		a->FillArcDensityHistogram(hist);
		//break;
	}
	
	return NOERROR;
}

//------------------------------------------------------------------
// FillSlopeIntDensityHistos
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FillSlopeIntDensityHistos(void)
{
	// Loop over track centers and find all hits which may contribute.
	// Remember all phi,z coordinates for these hits so we can histogram
	// slope and intercept of phi vs. z below.
	
	slope_density->Reset();
	offset_density->Reset();
	for(int i=0;i<Ndensity_histos;i++){
		slope_density_histos[i]->Reset();
		offset_density_histos[i]->Reset();
	}
	for(int j=0;j<Ncircles;j++){
		float phi[200], z[200];
		int Nhits = 0;
		DArcHit *a = archit;
		float x0 = circles[j].GetX1();
		float y0 = circles[j].GetY1();
		float r0 = sqrt(x0*x0 + y0*y0);
		float phi0 = atan2(y0, x0);
		if(phi0<0.0)phi0 += 2.0*M_PI;
		for(int i=0;i<Narchits;i++ ,a++){
			if(Nhits>=200)break;
			float d = a->DistToLine(x0,y0);
			if(d > masksize)continue;
			
			// Take cross-product to find relative angle between this hit
			// and vector pointing to origin
			float x = a->xhit+x0;
			float y = a->yhit+y0;
			float r = sqrt(x*x + y*y);
			float delta_phi = acos((x*x0 +y*y0)/(r*r0));
			
			phi[Nhits] = delta_phi;
			z[Nhits++] = a->zhit;
		}
		
		for(int i=0; i<Nhits-1; i++){
			for(int k=i+1; k<Nhits;k++){ 
				float m = (phi[i]-phi[k])/(z[i]-z[k]);
				float b = phi[i] - m*z[i];
				if(!finite(m) || !finite(b))continue;
				slope_density->Fill(m);
				offset_density->Fill(-b/m);
				if(j<Ndensity_histos){
					slope_density_histos[j]->Fill(m);
					offset_density_histos[j]->Fill(-b/m);
				}
			}
		}
	}
	
	return NOERROR;
}

//------------------------------------------------------------------
// DrawPhiZPoints
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::DrawPhiZPoints(void)
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
	
	// Since this is just for debugging, we can take the overhead of 
	// allocating and deleting every time we're called
	if(Nmarkers)delete markers;
	markers=NULL;
	Nmarkers = 0;
	markers = new TMarker[200];
	
	int colors[] = {kRed,kBlue,kMagenta,kGreen,kBlack};
	for(int j=0;j<Ncircles;j++){
		DArcHit *a = archit;
		float x0 = circles[j].GetX1();
		float y0 = circles[j].GetY1();
		float r0 = sqrt(x0*x0 + y0*y0);
		float phi0 = atan2(y0, x0);
		if(phi0<0.0)phi0 += 2.0*M_PI;
		for(int i=0;i<Narchits;i++ ,a++){
			
			float d = a->DistToLine(x0,y0);
			if(d > masksize)continue;
			
			// Take cross-product to find relative angle between this hit
			// and vector pointing to origin
			float x = a->xhit+x0;
			float y = a->yhit+y0;
			float r = sqrt(x*x + y*y);
			float delta_phi = acos((x*x0 +y*y0)/(r*r0));
			float z = a->zhit;
			
			markers[Nmarkers].SetX(z);
			markers[Nmarkers].SetY(delta_phi);
			markers[Nmarkers].SetMarkerColor(colors[j%6]);
			markers[Nmarkers].SetMarkerStyle(20+j);
			markers[Nmarkers].Draw();

			if(Nmarkers++>=200)break;
		}
		if(Nmarkers>=200)break;
	}

	return NOERROR;
}

//------------------------------------------------------------------
// SetNumDensityHistograms
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::SetNumDensityHistograms(int N)
{
	/// Set the number of density histograms to keep.
	///
	/// This is meant as a diagnostic and should not normally be used.
	/// It tells the factory to keep around multiple 2-D density
	/// histograms corresponding to subsequent stages of the track
	/// finding in the X/Y plane. Since the 2-D histogram can be
	/// quite large, setting this to a large number will consume
	/// a lot of memory (about 1/4 to 1/2 MB per histogram).
	/// The default value is set to 1 when the object is instantiated.
	/// The value of N here can be from 1 to 8.
	///
	/// Note also that this affects the number of slope_density
	/// and offset density histos also. Those are 1-D histos though
	/// so they do not take up us much space

	// Make sure N is in range
	if(N<1 || N>8){
		cerr<<__FILE__<<":"<<__LINE__<<" The number of density histograms"<<endl;
		cerr<<"requested is out of range. It must be between 1 and 8 inclusive."<<endl;
		cerr<<"The value passed was "<<N<<"."<<endl;
		
		return VALUE_OUT_OF_RANGE;
	}

	// If we're reducing the number of histos, delete the extras
	for(int i=N+1;i<=Ndensity_histos;i++){
		delete density_histos[i-1];
		delete slope_density_histos[i-1];
		delete offset_density_histos[i-1];
	}
	
	// If we're increasing the number of histos, instantiate them
	for(int i=Ndensity_histos;i<N;i++){
		density_histos[i] = new TH2F(*density);
		slope_density_histos[i] = new TH1F(*slope_density);
		offset_density_histos[i] = new TH1F(*offset_density);
	}
	
	Ndensity_histos = N;

	return NOERROR;
}

//------------------------------------------------------------------
// GetDensityHistogram
//------------------------------------------------------------------
TH2F* DFactory_MCTrackCandidates::GetDensityHistogram(int n)
{
	/// Return a pointer to 2-D density histogram n where n is value from
	/// 0 to 1 less than the last call to SetNumDensityHistograms().
	
	if(n<0 || n>=Ndensity_histos)return NULL;

	return density_histos[n];
}

//------------------------------------------------------------------
// GetSlopeDensityHistogram
//------------------------------------------------------------------
TH1F* DFactory_MCTrackCandidates::GetSlopeDensityHistogram(int n)
{
	/// Return a pointer to 1-D slope density histogram n where n is value from
	/// 0 to 1 less than the last call to SetNumDensityHistograms().
	
	if(n<0 || n>=Ndensity_histos)return NULL;

	return slope_density_histos[n];
}

//------------------------------------------------------------------
// GetOffsetDensityHistogram
//------------------------------------------------------------------
TH1F* DFactory_MCTrackCandidates::GetOffsetDensityHistogram(int n)
{
	/// Return a pointer to 1-D z-offset density histogram n where n is value from
	/// 0 to 1 less than the last call to SetNumDensityHistograms().
	
	if(n<0 || n>=Ndensity_histos)return NULL;

	return offset_density_histos[n];
}

//------------
// Print
//------------
derror_t DFactory_MCTrackCandidates::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	printheader("row:   Nhits:  x0(cm):   y0(cm):  z_vertex(cm):  dphi/dz(rad/cm):");
	
	MCTrackCandidate_t *trackcandidate = (MCTrackCandidate_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, trackcandidate++){

		printnewrow();
		
		printcol("%d",    i);
		printcol("%d",    trackcandidate->Nhits);
		printcol("%3.1f", trackcandidate->x0);
		printcol("%3.1f", trackcandidate->y0);
		printcol("%3.1f", trackcandidate->z_vertex);
		printcol("%1.3f", trackcandidate->dphidz);

		printrow();
	}
	cout<<endl;
}

