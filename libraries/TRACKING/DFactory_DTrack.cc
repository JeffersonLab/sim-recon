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
	
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		const DTrackCandidate *trackcandidate = trackcandidates[i];
		DTrack *track = new DTrack;
		
		double Ro = sqrt(trackcandidate->x0*trackcandidate->x0 + trackcandidate->y0*trackcandidate->y0);

		// The following is from a fit of Ro/sin(theta)/p_thrn vs. theta
		double par[] = {180.618178, -77.646922, 290.502314, -318.531368, 145.897184, -24.017713};
		double s = sin(trackcandidate->theta);
		double fun = par[0]+s*(par[1]+s*(par[2]+s*(par[3]+s*(par[4]+s*(par[5])))));

		track->q = trackcandidate->q;
		track->p = Ro/s/fun;
		track->theta = trackcandidate->theta;
		track->phi = trackcandidate->phi;
		track->x = 0.0;
		track->y = 0.0;
		track->z = trackcandidate->z_vertex;
		track->candidateid = trackcandidate->id;
		
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
		
		printcol("%x", i);
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
