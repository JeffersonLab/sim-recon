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

//------------
// Print
//------------
derror_t DFactory_FCALHits::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	cout<<name<<endl;
	cout<<"---------------------------------------"<<endl;
	cout<<"row:   x(cm):   y(cm):   E(GeV):   t(ns):"<<endl;
	cout<<endl;
	
	FCALHit_t *fcalhit = (FCALHit_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, fcalhit++){
		char str[80];
		memset(str,' ',80);
		str[79] = 0;

		char num[32];
		sprintf(num, "%d", i);
		strncpy(&str[3-strlen(num)], num, strlen(num));

		sprintf(num, "%3.1f", fcalhit->x);
		strncpy(&str[12-strlen(num)], num, strlen(num));
		sprintf(num, "%3.1f", fcalhit->y);
		strncpy(&str[21-strlen(num)], num, strlen(num));
		sprintf(num, "%2.3f", fcalhit->E);
		strncpy(&str[31-strlen(num)], num, strlen(num));
		sprintf(num, "%4.0f", fcalhit->t);
		strncpy(&str[40-strlen(num)], num, strlen(num));

		cout<<str<<endl;
	}
	cout<<endl;
}

