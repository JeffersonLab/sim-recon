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

//------------
// Print
//------------
derror_t DFactory_CDCHits::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	cout<<name<<endl;
	cout<<"---------------------------------------"<<endl;
	cout<<"row: radius(cm):  phim(rad):   dE(MeV):   t(ns):"<<endl;
	cout<<endl;
	
	CDCHit_t *cdchit = (CDCHit_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, cdchit++){
		char str[80];
		memset(str,' ',80);
		str[79] = 0;

		char num[32];
		sprintf(num, "%d", i);
		strncpy(&str[3-strlen(num)], num, strlen(num));

		sprintf(num, "%3.1f", cdchit->radius);
		strncpy(&str[15-strlen(num)], num, strlen(num));
		sprintf(num, "%1.3f", cdchit->phim);
		strncpy(&str[27-strlen(num)], num, strlen(num));
		sprintf(num, "%2.3f", cdchit->dE);
		strncpy(&str[38-strlen(num)], num, strlen(num));
		sprintf(num, "%4.0f", cdchit->t);
		strncpy(&str[47-strlen(num)], num, strlen(num));

		cout<<str<<endl;
	}
	cout<<endl;
}

