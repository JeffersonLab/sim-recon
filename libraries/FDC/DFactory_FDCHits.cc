// $Id$


#include "DFactory_FDCHits.h"


//-------------------
// evnt
//-------------------
derror_t DFactory_FDCHits::evnt(int eventnumber)
{

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(int i=0; i<PE->mult; i++){
		s_Chambers_t *chambers = NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardDC)
				chambers = PE->in[i].hitView->forwardDC->chambers;
		if(!chambers)continue;
		
		for(int j=0;j<chambers->mult;j++){
			int  layer = chambers->in[j].layer;
			int  module = chambers->in[j].module;

			s_AnodePlanes_t *anodeplanes = chambers->in[j].anodePlanes;
			if(anodeplanes){
			
				for(int k=0;k<anodeplanes->mult;k++){
					float tau = anodeplanes->in[k].tau;
					float z = anodeplanes->in[k].z;
					s_Wires_t *wires = anodeplanes->in[k].wires;
					if(!wires)continue;
				
					for(int m=0;m<wires->mult;m++){
						float u = wires->in[m].u;
						s_Hits_t *hits = wires->in[m].hits;
						if(!hits)continue;
						for(int n=0;n<hits->mult;n++){
							float dE = hits->in[n].dE;
							float t = hits->in[n].t;
					
							FDCHit_t *fdchit = (FDCHit_t*)_data->Add();
							fdchit->layer = layer;
							fdchit->module = module;
							fdchit->tau = tau;
							fdchit->z = z;
							fdchit->u = u;
							fdchit->dE = dE;
							fdchit->t = t;
							fdchit->type = 0;
						}
					}
				}
			}

			s_CathodePlanes_t *cathodeplanes = chambers->in[j].cathodePlanes;
			if(cathodeplanes){
			
				for(int k=0;k<cathodeplanes->mult;k++){
					float tau = cathodeplanes->in[k].tau;
					float z = cathodeplanes->in[k].z;
					s_Strips_t *strips = cathodeplanes->in[k].strips;
					if(!strips)continue;
				
					for(int m=0;m<strips->mult;m++){
						float u = strips->in[m].u;
					
						FDCHit_t *fdchit = (FDCHit_t*)_data->Add();
						fdchit->layer = layer;
						fdchit->module = module;
						fdchit->tau = tau;
						fdchit->z = z;
						fdchit->u = u;
						fdchit->dE = 0;
						fdchit->t = 0;
						fdchit->type = 1;
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
derror_t DFactory_FDCHits::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	cout<<name<<endl;
	cout<<"---------------------------------------"<<endl;
	cout<<"row: layer: module: tau(rad):    z(cm):  u(cm):  dE(MeV):   t(ns):   type:"<<endl;
	cout<<endl;
	
	FDCHit_t *fdchit = (FDCHit_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, fdchit++){
		char str[80];
		memset(str,' ',80);
		str[79] = 0;

		char num[32];
		sprintf(num, "%d", i);
		strncpy(&str[3-strlen(num)], num, strlen(num));

		sprintf(num, "%d", fdchit->layer);
		strncpy(&str[10-strlen(num)], num, strlen(num));
		sprintf(num, "%d", fdchit->module);
		strncpy(&str[18-strlen(num)], num, strlen(num));
		sprintf(num, "%3.1f", fdchit->tau);
		strncpy(&str[28-strlen(num)], num, strlen(num));
		sprintf(num, "%3.1f", fdchit->z);
		strncpy(&str[38-strlen(num)], num, strlen(num));
		sprintf(num, "%2.3f", fdchit->u);
		strncpy(&str[46-strlen(num)], num, strlen(num));
		sprintf(num, "%1.3f", fdchit->dE*1000.0);
		strncpy(&str[56-strlen(num)], num, strlen(num));
		sprintf(num, "%4.0f", fdchit->t);
		strncpy(&str[65-strlen(num)], num, strlen(num));
		sprintf(num, "%s", fdchit->type ? "anode":"cathode");
		strncpy(&str[73-strlen(num)], num, strlen(num));

		cout<<str<<endl;
	}
	cout<<endl;
}


