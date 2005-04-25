// $Id$
//
//    File: DFactory_DMCReconstructed.cc
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DMCTrackCandidate.h"
#include "DMCThrown.h"
#include "DFactory_DMCReconstructed.h"
#include "DEvent.h"

//------------------
// evnt
//------------------
derror_t DFactory_DMCReconstructed::evnt(int eventnumber)
{
	// For now, we just copy from the MCTrackCandidates. Eventually,
	// a track fitter will be implemented.
	vector<const DMCTrackCandidate*> mctc;
	event->Get(mctc);

	vector<const DMCThrown*> mcthrowns;
	event->Get(mcthrowns);
	
	for(int i=0; i<mctc.size(); i++){
		const DMCTrackCandidate *mctrackcandidate = mctc[i];
		DMCReconstructed *mcreconstructed = new DMCReconstructed;

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
		mcreconstructed->FindClosestThrown(mcthrowns);
		
		_data.push_back(mcreconstructed);
	}

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DMCReconstructed::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: type:  q:    p:    E: theta:   phi:   mass:     x:     y:     z:");
	
	for(int i=0; i<_data.size(); i++){

		DMCReconstructed *mcreconstructed = _data[i];

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
	
	return _table;
}
