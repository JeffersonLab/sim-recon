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
		s_HitView_t *HV = PE->in[i].hitView;
		s_Vcounters_t *vcounters = NULL;
		s_Hcounters_t *hcounters = NULL;
		if(PE->in[i].hitView){
			if(PE->in[i].hitView->forwardTOF){
				vcounters = PE->in[i].hitView->forwardTOF->vcounters;
				hcounters = PE->in[i].hitView->forwardTOF->hcounters;
			}
		}
		
		if(vcounters){
			for(int j=0;j<vcounters->mult;j++){
				float x = vcounters->in[j].x;
				s_Top_t *top = vcounters->in[j].top;
				if(top){
					s_Hits_t *hits = top->hits;
					for(int k=0;k<hits->mult;k++){
						float dE = hits->in[k].dE;
						float t = hits->in[k].t;
					
						TOFHit_t *tofhit = (TOFHit_t*)_data->Add();
						tofhit->x = x;
						tofhit->y = 0.0;
						tofhit->orientation = 0;
						tofhit->end = 0;
						tofhit->dE = dE;
						tofhit->t = t;
					}
				}
				s_Bottom_t *bottom = vcounters->in[j].bottom;
				if(bottom){
					s_Hits_t *hits = bottom->hits;
					for(int k=0;k<hits->mult;k++){
						float dE = hits->in[k].dE;
						float t = hits->in[k].t;
					
						TOFHit_t *tofhit = (TOFHit_t*)_data->Add();
						tofhit->x = x;
						tofhit->y = 0.0;
						tofhit->orientation = 0;
						tofhit->end = 1;
						tofhit->dE = dE;
						tofhit->t = t;
					}
				}
			}
		}

		if(hcounters){
			for(int j=0;j<hcounters->mult;j++){
				float y = hcounters->in[j].y;
				s_Left_t *left = hcounters->in[j].left;
				if(left){
					s_Hits_t *hits = left->hits;
					for(int k=0;k<hits->mult;k++){
						float dE = hits->in[k].dE;
						float t = hits->in[k].t;
					
						TOFHit_t *tofhit = (TOFHit_t*)_data->Add();
						tofhit->x = 0.0;
						tofhit->y = y;
						tofhit->orientation = 1;
						tofhit->end = 0;
						tofhit->dE = dE;
						tofhit->t = t;
					}
				}
				s_Right_t *right = hcounters->in[j].right;
				if(right){
					s_Hits_t *hits = right->hits;
					for(int k=0;k<hits->mult;k++){
						float dE = hits->in[k].dE;
						float t = hits->in[k].t;
					
						TOFHit_t *tofhit = (TOFHit_t*)_data->Add();
						tofhit->x = 0.0;
						tofhit->y = y;
						tofhit->orientation = 1;
						tofhit->end = 1;
						tofhit->dE = dE;
						tofhit->t = t;
					}
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

	printheader("row:   y(cm):      end:     dE(MeV):   t(ns):");
	
	TOFHit_t *tofhit = (TOFHit_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, tofhit++){

		printnewrow();
		
		printcol("%d", i);
		printcol("%3.1f", tofhit->y);
		printcol(tofhit->end ? "right":"left");
		printcol("%2.3fs", tofhit->dE*1000.0);
		printcol("%4.3fs", tofhit->t);
		
		printrow();
	}
	cout<<endl;

	return NOERROR;
}


