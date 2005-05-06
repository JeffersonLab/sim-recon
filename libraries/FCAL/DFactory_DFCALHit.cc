// $Id$
//
//    File: DFactory_DFCALHit.cc
// Created: Sun Apr  3 10:41:55 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DFactory_DFCALHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DFCALHit::evnt(int enventnumber)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_Rows_t *rows = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardEMcal)
				rows = PE->in[i].hitView->forwardEMcal->rows;
		if(!rows)continue;
		
		for(unsigned int j=0;j<rows->mult;j++){
			float y = rows->in[j].y;
			s_Columns_t *columns = rows->in[j].columns;
			if(!columns)continue;
			
			for(unsigned int k=0;k<columns->mult;k++){
				float x = columns->in[k].x;
				s_Showers_t *showers = columns->in[k].showers;
				if(!showers)continue;
				
				for(unsigned int m=0;m<showers->mult;m++){
					float E = showers->in[m].E;
					float t = showers->in[m].t;
					
					DFCALHit *fcalhit = new DFCALHit;
					fcalhit->x = x;
					fcalhit->y = y;
					fcalhit->E = E;
					fcalhit->t = t;
					_data.push_back(fcalhit);
				}
			}
		}
	}

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

	printheader("row:   x(cm):   y(cm):   E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALHit *fcalhit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%3.1f", fcalhit->x);
		printcol("%3.1f", fcalhit->y);
		printcol("%2.3f", fcalhit->E);
		printcol("%4.0f", fcalhit->t);
		printrow();
	}
	
	return _table;
}
