// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"
#include "hddm_s.h"


//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	eventNo = eventnumber;
	Ncdchits = 0;

	// Copy the cdc_trackhits into TVector3 array so
	// we can rotate the points.
	s_Cdc_trackhit_t *c = hddm->cdc_trackhits->in;
	TVector3 *v = cdchits;
	int *t = cdchit_tracks;
	for(Ncdchits=0;Ncdchits<hddm->cdc_trackhits->mult;Ncdchits++, c++, v++, t++){
		v->SetXYZ(c->x, c->y, c->z);
		*t = c->track;
	}


	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm->physicsEvents;
	if(!PE)return NOERROR;
	for(int i=0; i<PE->mult; i++){
		s_Reactions_t *reactions = PE->in[i].reactions;
		if(reactions){
			for(int j=0;j<reactions->mult;j++){
				s_Vertices_t *vertices = reactions->in[j].vertices;
				if(!vertices)continue;
				for(int k=0;k<vertices->mult;k++){
					s_Products_t *products = vertices->in[j].products;
					if(!products)continue;
					for(int m=0;m<products->mult;m++){
						if(products->in[m].momentum){
						}
					}
				}
			}
		}
	}

	return NOERROR;
}

