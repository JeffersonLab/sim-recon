// $Id$
//
//    File: DVertex_factory.cc
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <PID/DChargedTrack.h>

#include "DVertex_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DVertex_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DVertex_factory::brun(jana::JEventLoop *loop, int runnumber)
{
	// Initialize to reasonable defaults
	target_length = 0.0;
	target_z = 0.0;

	// Get Target parameters from XML
	DApplication *dapp = dynamic_cast<DApplication*> (loop->GetJApplication());
	DGeometry *geom = dapp ? dapp->GetDGeometry(runnumber):NULL;
	if(geom){
		geom->GetTargetLength(target_length);
		geom->GetTargetZ(target_z);
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DVertex_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get list of all charged tracks
	vector<const DChargedTrack*> chargedtracks;
	loop->Get(chargedtracks);
	
	// Copy most probable tracks to local container
	vector<const DTrackTimeBased*> trks;
	for(unsigned int i=0; i<chargedtracks.size(); i++){
		if(chargedtracks[i]->hypotheses.size()>0)trks.push_back(chargedtracks[i]->hypotheses[0]);
	}
	
	// Create a DVertex object and add it to our data list. The values
	// will be filled in below. This will need to be re-worked for detached
	// vertexes
	DVertex *vertex = new DVertex();
	_data.push_back(vertex);
	vertex->x.SetXYZ(0.0, 0.0, target_z);
	vertex->cov.ResizeTo(3,3);
	vertex->beamline_used = true;
	
	// If there are no tracks assume center of target which is already set so just return
	if(trks.size()==0)return NOERROR;
	
	// If there is only one track is found, use it's "position()" field which is the POCA to beamline
	if(trks.size()==1){
		vertex->x = trks[0]->position();
		vertex->AddAssociatedObject(trks[0]);
		return NOERROR;
	}
	
	// If we get here then there must be more than one track.
	
	// Simply average POCA to beamline for all tracks
	vertex->x.SetXYZ(0.0, 0.0, 0.0);
	for(unsigned int i=0; i<trks.size(); i++){
		const DTrackTimeBased *trk = trks[i];
		vertex->x += trk->position();
		vertex->AddAssociatedObject(trk);
	}
	vertex->x *= (1.0/(double)trks.size());
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DVertex_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DVertex_factory::fini(void)
{
	return NOERROR;
}

