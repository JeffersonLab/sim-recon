// $Id$


#include "DFactory_TOFHits.h"


//-------------------
// evnt
//-------------------
derror_t DFactory_TOFHits::evnt(int eventnumber)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(int i=0; i<PE->mult; i++){
		s_Slabs_t *slabs = NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardTOF)
				slabs = PE->in[i].hitView->forwardTOF->slabs;
		if(!slabs)continue;

		for(int j=0;j<slabs->mult;j++){
			float y = slabs->in[j].y;
			s_Left_t *left = slabs->in[j].left;
			if(left){
				s_Hits_t *hits = left->hits;
				if(!hits)continue;
				
				for(int m=0;m<hits->mult;m++){
					float dE = hits->in[m].dE;
					float t = hits->in[m].t;
					
					TOFHit_t *tofhit = (TOFHit_t*)_data->Add();
					tofhit->y = y;
					tofhit->end = 0;
					tofhit->dE = dE;
					tofhit->t = t;
				}
			}
			s_Right_t *right = slabs->in[j].right;
			if(right){
				s_Hits_t *hits = right->hits;
				if(!hits)continue;
				
				for(int m=0;m<hits->mult;m++){
					float dE = hits->in[m].dE;
					float t = hits->in[m].t;
					
					TOFHit_t *tofhit = (TOFHit_t*)_data->Add();
					tofhit->y = y;
					tofhit->end = 1;
					tofhit->dE = dE;
					tofhit->t = t;
				}
			}
		}
	}
	
	return NOERROR;

}

//-------------------
// Print
//-------------------
derror_t DFactory_TOFHits::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	cout<<name<<endl;
	cout<<"---------------------------------------"<<endl;
	cout<<"row:   y(cm):      end:     dE(MeV):   t(ns):"<<endl;
	cout<<endl;
	
	TOFHit_t *tofhit = (TOFHit_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, tofhit++){
		char str[80];
		memset(str,' ',80);
		str[79] = 0;

		char num[32];
		sprintf(num, "%d", i);
		strncpy(&str[3-strlen(num)], num, strlen(num));

		sprintf(num, "%3.1f", tofhit->y);
		strncpy(&str[12-strlen(num)], num, strlen(num));
		sprintf(num, "%s", tofhit->end ? "right":"left");
		strncpy(&str[22-strlen(num)], num, strlen(num));
		sprintf(num, "%2.3f", tofhit->dE*1000.0);
		strncpy(&str[35-strlen(num)], num, strlen(num));
		sprintf(num, "%4.3f", tofhit->t);
		strncpy(&str[44-strlen(num)], num, strlen(num));

		cout<<str<<endl;
	}
	cout<<endl;

	return NOERROR;
}


