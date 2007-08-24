// $Id: smear.cc 2432 2007-02-06 04:19:48Z davidl $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>
#include "HDDM/hddm_s.h"


//-----------
// Filter
//-----------
bool Filter(s_HDDM_t *hddm_s)
{
	// Return "true" to keep event, "false" to throw it away
	
	// Count number of charged and neutrals in final state
	int Ncharged = 0;
	int Nneutral = 0;
	
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return false;
	
	for(unsigned int i=0; i<PE->mult; i++){
		// ------------ Reactions --------------
		s_Reactions_t *reactions=PE->in[i].reactions;
		if(!reactions)continue;

		for(unsigned int j=0; j<reactions->mult; j++){
			s_Vertices_t *vertices = reactions->in[j].vertices;
			if(vertices){
				for(unsigned int k=0; k<vertices->mult; k++){
					s_Origin_t *origin = vertices->in[k].origin;
					s_Products_t *products = vertices->in[k].products;
					if(products && origin){
						for(unsigned int m=0;m<products->mult;m++){
							s_Product_t *product = &products->in[m];
							
							if(product->type<1)continue; // type==0 means intermediate state particle
							if(ParticleCharge(product->type)!=0.0){
								Ncharged++;
							}else{
								Nneutral++;
							}
						}
					}
				}
			}
		}
	}
	
	// Filter on there being at least 2 charged tracks
	if(Ncharged<2)return false;

	return true;
}
