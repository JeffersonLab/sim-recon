// $Id$


#include "DFactory_BCALHits.h"


//-------------------
// evnt
//-------------------
derror_t DFactory_BCALHits::evnt(int eventnumber)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(int i=0; i<PE->mult; i++){
		s_Modules_t *modules = NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->barrelEMcal)
				modules = PE->in[i].hitView->barrelEMcal->modules;
		if(!modules)continue;

		for(int j=0;j<modules->mult;j++){
			float phim = modules->in[j].phim;
			s_Upstream_t *upstream = modules->in[j].upstream;
			if(upstream){
				s_Showers_t *showers = upstream->showers;
				if(!showers)continue;
				
				for(int m=0;m<showers->mult;m++){
					float E = showers->in[m].E;
					float t = showers->in[m].t;
					
					BCALHit_t *bcalhit = (BCALHit_t*)_data->Add();
					bcalhit->phim = phim;
					bcalhit->end = 0;
					bcalhit->E = E;
					bcalhit->t = t;
				}
			}
			s_Downstream_t *downstream = modules->in[j].downstream;
			if(downstream){
				s_Showers_t *showers = downstream->showers;
				if(!showers)continue;
				
				for(int m=0;m<showers->mult;m++){
					float E = showers->in[m].E;
					float t = showers->in[m].t;
					
					BCALHit_t *bcalhit = (BCALHit_t*)_data->Add();
					bcalhit->phim = phim;
					bcalhit->end = 1;
					bcalhit->E = E;
					bcalhit->t = t;
				}
			}
		}
	}

	return NOERROR;
}

//-------------------
// Print
//-------------------
derror_t DFactory_BCALHits::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	cout<<name<<endl;
	cout<<"---------------------------------------"<<endl;
	cout<<"row:   phim(rad):      end:     E(GeV):   t(ns):"<<endl;
	cout<<endl;
	
	BCALHit_t *bcalhit = (BCALHit_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, bcalhit++){
		char str[80];
		memset(str,' ',80);
		str[79] = 0;

		char num[32];
		sprintf(num, "%d", i);
		strncpy(&str[3-strlen(num)], num, strlen(num));

		sprintf(num, "%1.3f", bcalhit->phim);
		strncpy(&str[16-strlen(num)], num, strlen(num));
		sprintf(num, "%s", bcalhit->end ? "downstream":"upstream");
		strncpy(&str[28-strlen(num)], num, strlen(num));
		sprintf(num, "%2.3f", bcalhit->E);
		strncpy(&str[38-strlen(num)], num, strlen(num));
		sprintf(num, "%4.0f", bcalhit->t);
		strncpy(&str[47-strlen(num)], num, strlen(num));

		cout<<str<<endl;
	}
	cout<<endl;

	return NOERROR;
}

