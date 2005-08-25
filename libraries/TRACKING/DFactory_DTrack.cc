// $Id$
//
//    File: DFactory_DTrack.cc
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DTrackCandidate.h"
#include "DFactory_DTrack.h"
#include "DEventLoop.h"

//------------------
// evnt
//------------------
derror_t DFactory_DTrack::evnt(DEventLoop *eventLoop, int eventnumber)
{
	// For now, we just copy from the MCTrackCandidates. Eventually,
	// a track fitter will be implemented.
	vector<const DTrackCandidate*> trackcandidates;
	eventLoop->Get(trackcandidates);
	
	identifier_t idcntr = 1;
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		const DTrackCandidate *trackcandidate = trackcandidates[i];
		DTrack *track = new DTrack;

		track->q = trackcandidate->q;
		track->p = trackcandidate->p;
		track->theta = trackcandidate->theta;
		track->phi = trackcandidate->phi;
		track->x = 0.0;
		track->y = 0.0;
		track->z = trackcandidate->z_vertex;
		track->candidateid = trackcandidate->id;
		track->id = idcntr++;
		
		_data.push_back(track);
	}

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DTrack::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: q:     p:   theta:   phi:    x:    y:    z:");
	
	for(unsigned int i=0; i<_data.size(); i++){

		DTrack *track = _data[i];

		printnewrow();
		
		printcol("%d", i);
		printcol("%+d", (int)track->q);
		printcol("%3.1f", track->p);
		printcol("%1.3f", track->theta);
		printcol("%1.3f", track->phi);
		printcol("%2.2f", track->x);
		printcol("%2.2f", track->y);
		printcol("%2.2f", track->z);

		printrow();
	}
	
	return _table;
}
