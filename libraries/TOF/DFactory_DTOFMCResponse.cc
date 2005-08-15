// $Id$
//
//    File: DFactory_DTOFMCResponse.cc
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#include "DFactory_DTOFMCResponse.h"

//------------------
// evnt
//------------------
derror_t DFactory_DTOFMCResponse::evnt(DEventLoop *loop, int eventnumber)
{
	// Code to generate factory data goes here. Add it like:
	//
	// DTOFMCResponse *myDTOFMCResponse = new DTOFMCResponse;
	// myDTOFMCResponse->x = x;
	// myDTOFMCResponse->y = y;
	// ...
	// _data.push_back(myDTOFMCResponse);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}


//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DTOFMCResponse::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
	v.clear();

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_Vcounters_t *vcounters = NULL;
		s_Hcounters_t *hcounters = NULL;
		if(PE->in[i].hitView){
			if(PE->in[i].hitView->forwardTOF){
				vcounters = PE->in[i].hitView->forwardTOF->vcounters;
				hcounters = PE->in[i].hitView->forwardTOF->hcounters;
			}
		}
		
		if(vcounters){
			for(unsigned int j=0;j<vcounters->mult;j++){
				float x = vcounters->in[j].x;
				float y = 0.0;
				s_Top_t *top = vcounters->in[j].top;
				if(top){
					s_Hits_t *hits = top->hits;
					for(unsigned int k=0;k<hits->mult;k++){
						float dE = hits->in[k].dE;
						float t = hits->in[k].t;
					
						DTOFMCResponse *tofhit = new DTOFMCResponse;
						tofhit->x = x;
						tofhit->y = y;
						tofhit->orientation = 0;
						tofhit->end = 0;
						tofhit->dE = dE;
						tofhit->t = t;
						v.push_back(tofhit);
					}
				}
				s_Bottom_t *bottom = vcounters->in[j].bottom;
				if(bottom){
					s_Hits_t *hits = bottom->hits;
					for(unsigned int k=0;k<hits->mult;k++){
						float dE = hits->in[k].dE;
						float t = hits->in[k].t;
					
						DTOFMCResponse *tofhit = new DTOFMCResponse;
						tofhit->x = x;
						tofhit->y = y;
						tofhit->orientation = 0;
						tofhit->end = 1;
						tofhit->dE = dE;
						tofhit->t = t;
						v.push_back(tofhit);
					}
				}
			}
		}

		if(hcounters){
			for(unsigned int j=0;j<hcounters->mult;j++){
				float y = hcounters->in[j].y;
				float x = 0.0;
				s_Left_t *left = hcounters->in[j].left;
				if(left){
					s_Hits_t *hits = left->hits;
					for(unsigned int k=0;k<hits->mult;k++){
						float dE = hits->in[k].dE;
						float t = hits->in[k].t;
					
						DTOFMCResponse *tofhit = new DTOFMCResponse;
						tofhit->x = x;
						tofhit->y = y;
						tofhit->orientation = 1;
						tofhit->end = 0;
						tofhit->dE = dE;
						tofhit->t = t;
						v.push_back(tofhit);
					}
				}
				s_Right_t *right = hcounters->in[j].right;
				if(right){
					s_Hits_t *hits = right->hits;
					for(unsigned int k=0;k<hits->mult;k++){
						float dE = hits->in[k].dE;
						float t = hits->in[k].t;
					
						DTOFMCResponse *tofhit = new DTOFMCResponse;
						tofhit->x = x;
						tofhit->y = y;
						tofhit->orientation = 1;
						tofhit->end = 1;
						tofhit->dE = dE;
						tofhit->t = t;
						v.push_back(tofhit);
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
const string DFactory_DTOFMCResponse::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	//		printheader("row:    x:     y:");
	//
	// 	for(int i=0; i<_data.size(); i++){
	//			DTOFMCResponse *myDTOFMCResponse = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDTOFMCResponse->x);
	//			printcol("%3.2f",	myDTOFMCResponse->y);
	//			printrow();
	//		}
	//
	return _table;

}
