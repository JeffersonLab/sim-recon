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
		Nqfit = 0;
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
	masksize2 = masksize*masksize;
	
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
	
	// Split tracks found in X/Y using the Z-info. This also fills in
	// the factory data
	FindLines();
	
	// At this point, we seem to often have tracks which were "found"
	// multiple times (i.e. they have essentially the same fitted
	// parameters). Someday, I'll have to track that down, but until
	// then, I'll just go through and weed out duplicates.
	MCTrackCandidate_t *a = (MCTrackCandidate_t*)_data->first();
	for(int i=0;i<_data->nrows-1; i++, a++){
		MCTrackCandidate_t *b = &a[1];
		int filtered = 0;
		for(int j=i+1;j<_data->nrows; j++, b++){
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
			a--;
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
// FindLines
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FindLines(void)
{
	/// Find lines in the phi-z plane and create tracks from them.
	///
	/// This gets called after FindCircles and uses the results
	/// of that to calculate phi about the circle centers for the
	/// hits belonging to the circle. A list of phi,z points is then
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

	int Nslope_density_bins = slope_density->GetXaxis()->GetNbins();
	
	// Delete last event's DQuickfit objects
	for(int i=0;i<Nqfit; i++)delete qfit[i];
	Nqfit = 0;

	// Loop over all circles
	for(int j=0;j<Ncircles;j++){

		// Fill the slope and z-intercept histos for all hits consistent
		// with a circle at x0,y0. At the same time, place the list of
		// hits in phi,z,Nhits.
		float x0 = circles[j].GetX1();
		float y0 = circles[j].GetY1();
		float r0 = sqrt(x0*x0 + y0*y0);
		FillSlopeIntDensityHistos(x0, y0);
		
		// For debugging
		if(j<Ndensity_histos){
			*slope_density_histos[j] = *slope_density;
			*offset_density_histos[j] = *offset_density;
		}

		// For now, assume a single vertex for all tracks included in
		// in this circle. Just use the maximum
		// in the z-intercept histo as the z coordinate of the vertex.
		int z_bin = offset_density->GetMaximumBin();
		float z_vertex = offset_density->GetBinCenter(z_bin);

		// Find all peaks in the slope histogram
		do{
			// Use a simple algorithm here. Just use the maximum as the
			// peak center. Zero out the histogram near the peak to look
			// for a new peak. The number of pairs of phi,z points
			// is N(N-1)/2 (I think?). Assuming at least 4 good hits are
			// needed for a track, we require that a maximum in the slope
			// histo be at least 6 units. In practice, 10 seems a better cut-off
			int slope_bin = slope_density->GetMaximumBin();
			if(slope_density->GetBinContent(slope_bin)<10)break;
			
			float slope = slope_density->GetBinCenter(slope_bin);
			
			// Create a new DQuickFit object to do the fit (should there
			// be enough hits to do one).
			DQuickFit *fit = new DQuickFit();
			
			// Here we add a row to our factory data since we want to fill
			// in the ihits values as we loop over the hits below. If it
			// turns out we don't have enough hits after looping over them
			// all, then we will delete this row.
			MCTrackCandidate_t *mctrackcandidate = (MCTrackCandidate_t*)_data->Add();
			mctrackcandidate->Nhits = 0;

			// At this point, the "used" flag in the array of DArcHits should
			// indicate which hits are consistent with the current circle
			// (i.e. circles[j]). We now look for all hits with the "used"
			// flag set that also have a delta_phi and zhit within a
			// certain distance of the line defined by slope and z_vertex.
			float m = slope;
			float b = -z_vertex*m;
			DArcHit *a = archit;
			for(int i=0;i<Narchits;i++ ,a++){
				if(!a->used)continue;

				// calculate distance to line squared.
				// NOTE: We assume here that we can always represent
				// phi as a function of z. This will NOT be the case for
				// tracks going out at 90 degrees from the target. This
				// will have to be fixed later.
				float z1 = (a->zhit-m*(b-a->delta_phi))/(1.0+m*m);
				float phi1 = m*z1 + b;
				float delta_z = a->zhit-z1;
				float delta_phi = r0*(a->delta_phi-phi1); // convert into cm
				float d2 = delta_z*delta_z + delta_phi*delta_phi;

				if(d2<masksize2){
					// Add hit to DQuickFit object
					fit->AddHit(a->rhit, a->phihit, a->zhit);
					
					// Add hit index to track in factory data
					mctrackcandidate->ihit[mctrackcandidate->Nhits++] = a->ihit;
				}
			}
			
			// If enough hits were added, then do the fit and record
			// the results
			if(fit->GetNhits()>=3){
				fit->FitTrack();
				mctrackcandidate->x0 = -fit->x0;	// why do we need the minus sign?
				mctrackcandidate->y0 = -fit->y0;	// why do we need the minus sign?
		
				mctrackcandidate->z_vertex = fit->z_vertex;
				mctrackcandidate->dphidz = fit->theta/r0;
				mctrackcandidate->p = fit->p;
				mctrackcandidate->p_trans = fit->p_trans;
				mctrackcandidate->q = fit->q;
				mctrackcandidate->phi = fit->phi;
				mctrackcandidate->theta = fit->theta;
			}else{
				// Oops! not enough hits for a track. Delete this one.
				_data->Delete(_data->nrows-1);
			}
			
			// Keep the DQuickFit object around (unless we have too many)
			if(Nqfit<32){
				qfit[Nqfit++] = fit;
			}else{
				delete fit;
			}

			// Zero out the slope_density histo for 5 bins on either side
			// of the peak.
			for(int i=slope_bin-5;i<=slope_bin+5;i++){
				if(i<1 || i>Nslope_density_bins)continue;
				slope_density->SetBinContent(i, 0.0);
			}

		}while(1);

	}
	
	// For debugging
	slope_density->Reset();
	offset_density->Reset();
	for(int i=0; i<Ndensity_histos; i++){
		slope_density->Add(slope_density_histos[i]);
		offset_density->Add(offset_density_histos[i]);
	}

	return NOERROR;
}

//------------------------------------------------------------------
// FillSlopeIntDensityHistos
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FillSlopeIntDensityHistos(float x0, float y0)
{
	/// This is called by FindLines()
	
	slope_density->Reset();
	offset_density->Reset();
	float r0 = sqrt(x0*x0 + y0*y0);
	float x0_unit = x0/r0;
	float y0_unit = y0/r0;
	
	// Loop over all hits. Use the "used" flag in the DArchit objects to
	// record who is and isn't used to fill the slope and z-intercept histos
	DArcHit *a = archit;
	for(int i=0;i<Narchits;i++ ,a++){
		float d2 = a->Dist2ToLine(x0,y0);
		if(d2 > masksize2){
			a->used = 0;
			continue;
		}
		a->used = 1;

		// Take cross-product to find relative angle between this hit
		// and vector pointing to origin
		float x = a->xhit+x0;
		float y = a->yhit+y0;
		float r = sqrt(x*x + y*y);
		a->delta_phi = acos((x*x0_unit +y*y0_unit)/r);
	}

	// Loop over all pairs of phi,z points, filling the slope
	// and z-intercept histos
	a = archit;
	for(int i=0; i<Narchits-1; i++, a++){
		if(!a->used)continue;
		DArcHit *b = &a[1];
		for(int k=i+1; k<Narchits; k++, b++){ 
			if(!b->used)continue;
			float m = (a->delta_phi - b->delta_phi)/(a->zhit - b->zhit);
			float z = a->delta_phi - m*a->zhit;
			if(!finite(m) || !finite(z))continue;
			slope_density->Fill(m);
			offset_density->Fill(-z/m);
		}
	}
	
	return NOERROR;
}

//------------------------------------------------------------------
// GetQFit
//------------------------------------------------------------------
DQuickFit* DFactory_MCTrackCandidates::GetQFit(int n)
{
	if(n<0 || n>=Nqfit)return NULL;
	
	return qfit[n];
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

//------------------------------------------------------------------
// ThereCanBeOnlyOne
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::ThereCanBeOnlyOne(int trk1, int trk2)
{
	/// See the comment at the end of evnt(). Basically, we need to choose
	/// which one track to keep and adjust the _data array to keep it.
	
	// first, make sure trk1 and trk2 are in sequential order
	int trka = trk1<trk2 ? trk1:trk2;
	int trkb = trk1<trk2 ? trk2:trk1;

	// If we're keeping the track further down the list, then copy
	// its results into the position near the front of the list.
	// We should be using the chisq from the fits, but that isn't
	// being calculated for the 3D track at the moment.
	if(qfit[trka]->GetNhits() < qfit[trkb]->GetNhits()){
		MCTrackCandidate_t *a = (MCTrackCandidate_t*)_data->index(trka);
		MCTrackCandidate_t *b = (MCTrackCandidate_t*)_data->index(trkb);
		*a = *b;
		DQuickFit *tmp = qfit[trka];
		qfit[trka] = qfit[trkb];
		qfit[trkb] = tmp;
	}
	
	// Shift the tracks after trkb up in the list(s)
	_data->Delete(trkb);
	delete qfit[trkb];
	Nqfit--;
	for(int i=trkb; i<Nqfit; i++)qfit[i] = qfit[i+1];

	return NOERROR;	
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

	printheader("row: Nhits: x0(cm): y0(cm): z_vertex: dphi/dz:  q:   p: p_trans:   phi: theta:");
	
	MCTrackCandidate_t *trackcandidate = (MCTrackCandidate_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, trackcandidate++){

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
	cout<<endl;
}



//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// ---------- Things below are unused -----------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


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
