// $Id$


#include "DFactory_CDCHits.h"


//-------------------
// evnt
//-------------------
derror_t DFactory_CDCHits::evnt(int eventnumber)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(int i=0; i<PE->mult; i++){
		// ------------ CdcPoints, Hits --------------
		s_Rings_t *rings=NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->centralDC)
				rings = PE->in[i].hitView->centralDC->rings;
		if(rings){
			for(int j=0;j<rings->mult;j++){
				float radius = rings->in[j].radius;
				s_Straws_t *straws = rings->in[j].straws;
				if(straws){
					for(int k=0;k<straws->mult;k++){
						float phim = straws->in[k].phim;
						s_Hits_t *hits = straws->in[k].hits;
						if(hits){
							for(int m=0;m<hits->mult;m++){
								float dE = hits->in[m].dE;
								float t = hits->in[m].t;
								
								// Add a row to the factory data
								CDCHit_t *cdchit = (CDCHit_t*)_data->Add();
								cdchit->radius = radius;
								cdchit->phim = phim;
								cdchit->dE = dE;
								cdchit->t = t;
							}
						}
					}
				}
			}
		}
	}

	
	return NOERROR;
}


