// $Id$

#include "DEvent.h"
#include "DFactory_MCThrown.h"
#include "particleType.h"

//-------------------
// evnt
//-------------------
derror_t DFactory_MCThrown::evnt(int eventnumber)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(int i=0; i<PE->mult; i++){
		// ------------ Reactions --------------
		s_Reactions_t *reactions=PE->in[i].reactions;
		if(!reactions)continue;

		for(int j=0; j<reactions->mult; j++){
			s_Vertices_t *vertices = reactions->in[j].vertices;
			if(vertices){
				for(int k=0; k<vertices->mult; k++){
					s_Origin_t *origin = vertices->in[k].origin;
					s_Products_t *products = vertices->in[k].products;
					if(products && origin){
						for(int m=0;m<products->mult;m++){
							s_Product_t *product = &products->in[m];
							
							MCThrown_t *mcthrown = (MCThrown_t*)_data->Add();
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
						}
					}
				}
			}
		}
	}

	return NOERROR;
}

//------------
// Print
//------------
derror_t DFactory_MCThrown::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	printheader("row: type:  q:    p:    E: theta:   phi:   mass:     x:     y:     z:");
	
	MCThrown_t *mcthrown = (MCThrown_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, mcthrown++){

		printnewrow();
		
		printcol("%d", i);
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
	cout<<endl;
}

