// $Id$
//
//    File: DBCALHit_factory.cc
// Created: Thu Jun  9 10:14:35 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DBCALHit_factory.h"

//------------------
// evnt
//------------------
jerror_t DBCALHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	/// Place holder for now. 

	return NOERROR;
}

//------------------
// Extract_HDDM
//------------------
jerror_t DBCALHit_factory::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
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
			hits->barrelEMcal == HDDM_NULL ||
			hits->barrelEMcal->bcalCells == HDDM_NULL)continue;
		
		// Loop over BCAL cells
		s_BcalCells_t *cells = hits->barrelEMcal->bcalCells;
		for(unsigned int j=0;j<cells->mult;j++){
			s_BcalCell_t *cell = &cells->in[j];
			if(cell->bcalUpstreamHits != HDDM_NULL){
				for(unsigned int k=0; k<cell->bcalUpstreamHits->mult; k++){
					s_BcalUpstreamHit_t *upstreamhit = &cell->bcalUpstreamHits->in[k];
					
					DBCALHit *bcalhit = new DBCALHit();
					bcalhit->module = cell->module;
					bcalhit->layer = cell->layer;
					bcalhit->sector = cell->sector;
					bcalhit->end = DBCALHit::UPSTREAM;
					bcalhit->E = upstreamhit->E;
					bcalhit->t = upstreamhit->t;
					v.push_back(bcalhit);
				}
			}

			if(cell->bcalDownstreamHits != HDDM_NULL){
				for(unsigned int k=0; k<cell->bcalDownstreamHits->mult; k++){
					s_BcalDownstreamHit_t *downstreamhit = &cell->bcalDownstreamHits->in[k];
					
					DBCALHit *bcalhit = new DBCALHit();
					bcalhit->module = cell->module;
					bcalhit->layer = cell->layer;
					bcalhit->sector = cell->sector;
					bcalhit->end = DBCALHit::DOWNSTREAM;
					bcalhit->E = downstreamhit->E;
					bcalhit->t = downstreamhit->t;
					v.push_back(bcalhit);
				}
			}
		} // j   (cells)
	} // i   (physicsEvents)

	return NOERROR;
}

//------------------
// toString
//------------------
const string DBCALHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()==0)return string(); // don't print anything if we have no data!

	printheader("row:   module:  layer:  sector:         end:     E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DBCALHit *bcalhit = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%d",	bcalhit->module);
		printcol("%d",	bcalhit->layer);
		printcol("%d",	bcalhit->sector);
		printcol(bcalhit->end==DBCALHit::UPSTREAM ? "upstream":"downstream");
		printcol("%2.3f",	bcalhit->E);
		printcol("%4.0f",	bcalhit->t);
		printrow();
	}

	return _table;
}
