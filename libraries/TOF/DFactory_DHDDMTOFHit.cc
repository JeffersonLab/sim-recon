// $Id$
//
//    File: DFactory_DHDDMTOFHit.cc
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>	

#include "DFactory_DHDDMTOFHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DHDDMTOFHit::evnt(DEventLoop *loop, int eventnumber)
{
	// no code should be here -- this factory is used strictly for reading in
	// HDDM data
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DHDDMTOFHit::toString(void)
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
	//			DHDDMTOFHit *myDHDDMTOFHit = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDHDDMTOFHit->x);
	//			printcol("%3.2f",	myDHDDMTOFHit->y);
	//			printrow();
	//		}
	//
	return _table;

}




//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DHDDMTOFHit::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
    v.clear();

	// Loop over Physics Events
    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE) return NOERROR;
	
    identifier_t identifier = 0;

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

                s_Top_t *top = vcounters->in[j].top;
                s_Bottom_t *bottom = vcounters->in[j].bottom;

                float y = vcounters->in[j].x;

                if(top){
                    s_Hits_t *hits = top->hits;
                    for(unsigned int k=0;k<hits->mult;k++){		
                        DHDDMTOFHit *tofhit = new DHDDMTOFHit;
                        tofhit->orientation = 0;
                        tofhit->end         = 0;
                        tofhit->y           = y;
                        tofhit->t           = hits->in[k].t;
                        tofhit->E           = hits->in[k].dE;
                        identifier++;
                        v.push_back(tofhit);
                    }
                }

                if(bottom){
                    s_Hits_t *hits = bottom->hits;
                    for(unsigned int k=0;k<hits->mult;k++){		
                        DHDDMTOFHit *tofhit = new DHDDMTOFHit;
                        tofhit->orientation = 0;
                        tofhit->end         = 1;
                        tofhit->y           = y;
                        tofhit->t           = hits->in[k].t;
                        tofhit->E           = hits->in[k].dE;
                        identifier++;
                        v.push_back(tofhit);
                    }
                }

            }
        }


        if(hcounters){

            for(unsigned int j=0;j<hcounters->mult;j++){

                s_Left_t *left = hcounters->in[j].left;
                s_Right_t *right = hcounters->in[j].right;

                float y = hcounters->in[j].y;

                if(left){
                    s_Hits_t *hits = left->hits;
                    for(unsigned int k=0;k<hits->mult;k++){		
                        DHDDMTOFHit *tofhit = new DHDDMTOFHit;
                        tofhit->orientation = 1;
                        tofhit->end         = 0;
                        tofhit->y           = y;
                        tofhit->t           = hits->in[k].t;
                        tofhit->E           = hits->in[k].dE;
                        identifier++;
                        v.push_back(tofhit);
                    }
                }

                if(right){
                    s_Hits_t *hits = right->hits;
                    for(unsigned int k=0;k<hits->mult;k++){		
                        DHDDMTOFHit *tofhit = new DHDDMTOFHit;
                        tofhit->orientation = 1;
                        tofhit->end         = 1;
                        tofhit->y           = y;
                        tofhit->t           = hits->in[k].t;
                        tofhit->E           = hits->in[k].dE;
                        identifier++;
                        v.push_back(tofhit);
                    }
                }

            }
        }

    }

    return NOERROR;

}


