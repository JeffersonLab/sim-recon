// $Id$
//
//    File: DFactory_DFCALHit.cc
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <cassert>

#include "DFactory_DFCALHit.h"
#include "DFCALHit.h"
#include "DFCALMCResponse.h"
#include "DFCALGeometry.h"


//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DFCALHit::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	// extract the FCAL Geometry
	vector<const DFCALGeometry*> fcalGeomVect;
	eventLoop->Get( fcalGeomVect );
	if(fcalGeomVect.size() != 1){
		cerr<<__FILE__<<":"<<__LINE__<<" Bad number of DFCALGeometry objects ("<<fcalGeomVect.size()<<")!"<<endl;
		return VALUE_OUT_OF_RANGE;
	}
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardEMcal == HDDM_NULL ||
			hits->forwardEMcal->fcalBlocks == HDDM_NULL)continue;

		s_FcalBlocks_t *blocks = hits->forwardEMcal->fcalBlocks;
		for(unsigned int j=0; j<blocks->mult; j++){
			s_FcalBlock_t *block = &blocks->in[j];
			for(unsigned int k=0; k<block->fcalHits->mult; k++){
				s_FcalHit_t *fcalhit = &block->fcalHits->in[k];
				TVector2 position = fcalGeom.positionOnFace(block->row, block->column);
				
				DFCALHit *dfcalhit = new DFCALHit();
				dfcalhit->column = block->column;
				dfcalhit->row = block->row;
				dfcalhit->x = position.X();
				dfcalhit->y = position.Y();
				dfcalhit->E = fcalhit->E;
				dfcalhit->t = fcalhit->t;
				
				v.push_back((void*)dfcalhit);

			} // k  (fcalhits)
		} // j  (blocks)
	} // i  (physicsEvents)
	
	return NOERROR;
}


//------------------
// toString
//------------------
const string DFactory_DFCALHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:  column:   row:   x(cm):   y(cm):   E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALHit *fcalhit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%d", fcalhit->column);
		printcol("%d", fcalhit->row);
		printcol("%3.1f", fcalhit->x);
		printcol("%3.1f", fcalhit->y);
		printcol("%2.3f", fcalhit->E);
		printcol("%4.0f", fcalhit->t);
		printrow();
	}
	
	return _table;
}
