// $Id$
//
//    File: DFactory_DCDCHit.cc
// Created: Thu Jun  9 10:22:37 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DFactory_DCDCHit.h"

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DCDCHit::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
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
			hits->centralDC == HDDM_NULL ||
			hits->centralDC->cdcStraws == HDDM_NULL)continue;

		s_CdcStraws_t *straws = hits->centralDC->cdcStraws;
		for(unsigned int j=0; j<straws->mult; j++){
			s_CdcStraw_t *straw = &straws->in[j];
			for(unsigned int k=0; k<straw->cdcStrawHits->mult; k++){
				s_CdcStrawHit_t *strawHit = &straw->cdcStrawHits->in[k];
				
				// Add a row to the factory data
				DCDCHit *cdchit = new DCDCHit;
				cdchit->ring = straw->ring;
				cdchit->straw = straw->straw;
				cdchit->radius = 0.0; // need to get this from ring,straw
				cdchit->phim = 0.0; // need to get this from ring,straw
				cdchit->dE = strawHit->dE;
				cdchit->t = strawHit->t;
				v.push_back(cdchit);
			} // k (strawHits)
		} // j (straws)
	} // i  (physicsEvents)

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DCDCHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: ring:  straw:  radius(cm):  phim(rad):   dE(MeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DCDCHit *cdchit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%d", cdchit->ring);
		printcol("%d", cdchit->straw);
		printcol("%3.1f", cdchit->radius);
		printcol("%1.3f", cdchit->phim);
		printcol("%2.3f", cdchit->dE);
		printcol("%4.0f", cdchit->t);
		printrow();
	}

	return _table;
}
