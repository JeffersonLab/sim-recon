// $Id$
//
//    File: DHDDMTOFHit_factory.cc
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>	

#include "DHDDMTOFHit_factory.h"

//------------------
// evnt
//------------------
jerror_t DHDDMTOFHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// no code should be here -- this factory is used strictly for reading in
	// HDDM data
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DHDDMTOFHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   paddle: plane: end:    t:        dE:      id:");

	for(unsigned int i=0; i<_data.size(); i++){
		DHDDMTOFHit *hit = _data[i];

		printnewrow();
		printcol("%d",		hit->paddle);
		printcol("%d",		hit->plane);
		printcol("%d",		hit->end);
		printcol("%1.3f",	hit->t);
		printcol("%1.3f",	hit->dE);
		printcol("%x",		hit->id);
		printrow();
	}

	return _table;

}




//------------------
// Extract_HDDM
//------------------
jerror_t DHDDMTOFHit_factory::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects.
	
    v.clear();

	// Loop over Physics Events
    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardTOF == HDDM_NULL ||
			hits->forwardTOF->ftofCounters == HDDM_NULL)continue;
		
		s_FtofCounters_t* ftofCounters = hits->forwardTOF->ftofCounters;
		
		// Loop over counters
		s_FtofCounter_t *ftofCounter = ftofCounters->in;
		for(unsigned int j=0;j<ftofCounters->mult; j++, ftofCounter++){
			 
			// Loop over left hits
			s_FtofLeftHits_t *ftofLeftHits = ftofCounter->ftofLeftHits;
			s_FtofLeftHit_t *ftofLeftHit = ftofLeftHits->in;
			for(unsigned int k=0;k<ftofLeftHits->mult; k++, ftofLeftHit++){
				DHDDMTOFHit *tofhit = new DHDDMTOFHit;
				tofhit->paddle	= ftofCounter->paddle;
				tofhit->plane	= ftofCounter->plane;
				tofhit->end		= 0;
				tofhit->t		= ftofLeftHit->t;
				tofhit->dE		= ftofLeftHit->dE;
				v.push_back(tofhit);
			}
			 
			 // Loop over right hits
			s_FtofRightHits_t *ftofRightHits = ftofCounter->ftofRightHits;
			s_FtofRightHit_t *ftofRightHit = ftofRightHits->in;
			for(unsigned int k=0;k<ftofRightHits->mult; k++, ftofRightHit++){
				DHDDMTOFHit *tofhit = new DHDDMTOFHit;
				tofhit->paddle	= ftofCounter->paddle;
				tofhit->plane	= ftofCounter->plane;
				tofhit->end		= 1;
				tofhit->t		= ftofRightHit->t;
				tofhit->dE		= ftofRightHit->dE;
				v.push_back(tofhit);
			}
		}
	}

    return NOERROR;

}


