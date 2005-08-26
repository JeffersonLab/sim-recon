// $Id$
//
//    File: DFactory_DMCThrown.cc
// Created: Sun Apr  3 12:22:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DFactory_DMCThrown.h"

//------------------
// evnt
//------------------
derror_t DFactory_DMCThrown::evnt(DEventLoop *loop, int enventnumber)
{
	/// This doesn't do anything. All of the work is done in  Extract_HDDM()

	return NOERROR;
}

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DMCThrown::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
	v.clear();
	identifier_t idcntr = 1;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
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
							
							DMCThrown *mcthrown = new DMCThrown;
							mcthrown->x = origin->vx;
							mcthrown->y = origin->vy;
							mcthrown->z = origin->vz;
							mcthrown->type = product->type;
							mcthrown->q = (float)ParticleCharge(product->type);
							mcthrown->mass = ParticleMass(product->type);
							mcthrown->E = product->momentum->E;
							mcthrown->p = sqrt(mcthrown->E*mcthrown->E - mcthrown->mass*mcthrown->mass);
							mcthrown->phi = atan2(product->momentum->py, product->momentum->px);
							if(mcthrown->phi<0.0)mcthrown->phi += 2.0*M_PI;
							mcthrown->theta = acos(product->momentum->pz/mcthrown->p);
							mcthrown->id = idcntr++;
							v.push_back((void*)mcthrown);
						}
					}
				}
			}
		}
	}
	

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DMCThrown::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("id: type:  q:    p:    E: theta:   phi:   mass:     x:     y:     z:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DMCThrown * mcthrown = _data[i];

		printnewrow();
		
		printcol("%d", mcthrown->id);
		printcol("%d", mcthrown->type);
		printcol("%+d", (int)mcthrown->q);
		printcol("%3.1f", mcthrown->p);
		printcol("%3.1f", mcthrown->E);
		printcol("%1.3f", mcthrown->theta);
		printcol("%1.3f", mcthrown->phi);
		printcol("%1.3f", mcthrown->mass);
		printcol("%2.2f", mcthrown->x);
		printcol("%2.2f", mcthrown->y);
		printcol("%2.2f", mcthrown->z);

		printrow();
	}
	
	return _table;
}
