// $Id$

#include "DEvent.h"
#include "DFactory_MCReconstructed.h"
#include "DFactory_MCTrackCandidates.h"
#include "particleType.h"

//-------------------
// evnt
//-------------------
derror_t DFactory_MCReconstructed::evnt(int eventnumber)
{
	// For now, we just copy from the MCTrackCandidates. Eventually,
	// a track fitter will be implemented.
	DContainer *mctc = event->Get("MCTrackCandidates");
	
	MCTrackCandidate_t *mctrackcandidate = (MCTrackCandidate_t*)mctc->first();
	for(int i=0; i<mctc->nrows; i++, mctrackcandidate++){
		MCReconstructed_t *mcreconstructed = (MCReconstructed_t*)_data->Add();

		mcreconstructed->type = 0;
		mcreconstructed->q = mctrackcandidate->q;
		mcreconstructed->p = mctrackcandidate->p;
		mcreconstructed->E = 0.0;
		mcreconstructed->theta = mctrackcandidate->theta;
		mcreconstructed->phi = mctrackcandidate->phi;
		mcreconstructed->x = 0.0;
		mcreconstructed->y = 0.0;
		mcreconstructed->z = mctrackcandidate->z_vertex;
		mcreconstructed->mass = 0.0;
	}

	return NOERROR;
}

//------------
// Print
//------------
derror_t DFactory_MCReconstructed::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	printheader("row: type:  q:    p:    E: theta:   phi:   mass:     x:     y:     z:");
	
	MCReconstructed_t *mcreconstructed = (MCReconstructed_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, mcreconstructed++){

		printnewrow();
		
		printcol("%d", i);
		printcol("%d", mcreconstructed->type);
		printcol("%+d", (int)mcreconstructed->q);
		printcol("%3.1f", mcreconstructed->p);
		printcol("%3.1f", mcreconstructed->E);
		printcol("%1.3f", mcreconstructed->theta);
		printcol("%1.3f", mcreconstructed->phi);
		printcol("%1.3f", mcreconstructed->mass);
		printcol("%2.2f", mcreconstructed->x);
		printcol("%2.2f", mcreconstructed->y);
		printcol("%2.2f", mcreconstructed->z);

		printrow();
	}
	cout<<endl;
}

