// $Id$


#include "DFactory_FCALHits.h"


//-------------------
// evnt
//-------------------
derror_t DFactory_FCALHits::evnt(int eventnumber)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(int i=0; i<PE->mult; i++){
		s_Rows_t *rows = NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardEMcal)
				rows = PE->in[i].hitView->forwardEMcal->rows;
		if(!rows)continue;
		
		for(int j=0;j<rows->mult;j++){
			float y = rows->in[j].y;
			s_Columns_t *columns = rows->in[j].columns;
			if(!columns)continue;
			
			for(int k=0;k<columns->mult;k++){
				float x = columns->in[k].x;
				s_Showers_t *showers = columns->in[k].showers;
				if(!showers)continue;
				
				for(int m=0;m<showers->mult;m++){
					float E = showers->in[m].E;
					float t = showers->in[m].t;
					
					FCALHit_t *fcalhit = (FCALHit_t*)_data->Add();
					fcalhit->x = x;
					fcalhit->y = y;
					fcalhit->E = E;
					fcalhit->t = t;
				}
			}
		}
	}
	
	return NOERROR;
}


