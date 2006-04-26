// $Id$
//
//    File: DFactory_DHDDMTOFTruth.cc
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>	

#include "DFactory_DHDDMTOFTruth.h"

//------------------
// evnt
//------------------
derror_t DFactory_DHDDMTOFTruth::evnt(DEventLoop *loop, int eventnumber)
{
	// no code should be here -- this factory is used strictly for reading in
	// HDDM data
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DHDDMTOFTruth::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:  primary: track: t:       x:     y:     z:      orientation:");

	for(unsigned int i=0; i<_data.size(); i++){
		DHDDMTOFTruth *truth = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%d",	truth->primary);
		printcol("%d",	truth->track);
		printcol("%1.3f",	truth->t);
		printcol("%1.3f",	truth->x);
		printcol("%1.3f",	truth->y);
		printcol("%1.3f",	truth->z);
		printcol("%d",	truth->orientation);
		printrow();
	}

	return _table;

}




//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DHDDMTOFTruth::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{

	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
    v.clear();
	
        // Loop over Physics Events
    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE) return NOERROR;
	
    for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardTOF == HDDM_NULL ||
			hits->forwardTOF->ftofTruthPoints == HDDM_NULL)continue;

		s_FtofTruthPoints_t* ftofTruthPoints = hits->forwardTOF->ftofTruthPoints;

		// Loop truth hits
		s_FtofTruthPoint_t *ftofTruthPoint = ftofTruthPoints->in;
		for(unsigned int j=0;j<ftofTruthPoints->mult; j++, ftofTruthPoint++){
			DHDDMTOFTruth *toftruth = new DHDDMTOFTruth;
		
			toftruth->orientation = 1;
			toftruth->primary     = ftofTruthPoint->primary;
			toftruth->t           = ftofTruthPoint->t;
			toftruth->track       = ftofTruthPoint->track;
			toftruth->x           = ftofTruthPoint->x;
			toftruth->y           = ftofTruthPoint->y;
			toftruth->z           = ftofTruthPoint->z;
			v.push_back(toftruth);
		}
	}

	return NOERROR;
}

