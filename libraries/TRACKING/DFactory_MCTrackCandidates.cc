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
	
	// set limits for plot. This represents the space where the center 
	// of the circle can be. It can be (and often is) outside of the
	// bounds of the solenoid.
	circle_max = 120.0; // in cm.
	
	// The number of bins per cm (in one dimension) for the density histogram
	bins_per_cm = 1.0;

	int Nbins = (int)round(2.0*circle_max/bins_per_cm);
	density = new TH2F("density","Density",Nbins,-circle_max,circle_max,Nbins,-circle_max,circle_max);	
	Ndensity_histos = 1;
	density_histos[0] = density;

	// max distance a line-of-circle-centers line can
	// be from a focal point and still be considered on the circle
	masksize = 5.0; // in cm
	
	// See note in FindCircles
	flip_x_axis = 0;
}

//------------------------------------------------------------------
// ~DFactory_MCTrackCandidates (destructor)
//------------------------------------------------------------------
DFactory_MCTrackCandidates::~DFactory_MCTrackCandidates()
{
	for(int i=0; i<Ndensity_histos; i++)delete density_histos[i];
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
	for(Narchits=0;Narchits<mccheathits->nrows;Narchits++, mccheathit++){
	
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
		
		archit[Narchits].SetXYZ(x,y,mccheathit->z);
	}
	
	// Find circle patterns first. (The results are left in circles[])
	FindCircles();
		
	return NOERROR;
}

//------------------------------------------------------------------
// FindTracks
//------------------------------------------------------------------
derror_t DFactory_MCTrackCandidates::FindCircles(void)
{
	/// Find circles by repeatedly filling the density
	/// histogram and looking for the maximum. As each
	/// maximum is found, record its position in the
	/// circles[] array.

	Ncircles = 0;	
	do{
		// use tmp so we can more easily build in a mechanism for
		// filling multiple histograms
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
			if(d <= masksize){
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

	// Make sure N is in range
	if(N<1 || N>8){
		cerr<<__FILE__<<":"<<__LINE__<<" The number of density histograms"<<endl;
		cerr<<"requested is out of range. It must be between 1 and 8 inclusive."<<endl;
		cerr<<"The value passed was "<<N<<"."<<endl;
		
		return VALUE_OUT_OF_RANGE;
	}

	// If we're reducing the number of histos, delete the extras
	for(int i=N+1;i<=Ndensity_histos;i++)delete density_histos[i-1];
	
	// If we're increasing the number of histos, instantiate them
	for(int i=Ndensity_histos;i<N;i++)density_histos[i] = new TH2F(*density);
	
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

