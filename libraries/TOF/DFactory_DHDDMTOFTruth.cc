// $Id$
//
//    File: DFactory_DHDDMTOFTruth.cc
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>	

#include "DFactory_DHDDMTOFTruth.h"

//------------------
// evnt
//------------------
derror_t DFactory_DHDDMTOFTruth::evnt(DEventLoop *loop, int eventnumber)
{
	// no code should be here -- this factory is used strictly for reading in
	// HDDM data
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DHDDMTOFTruth::toString(void)
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
	//			DHDDMTOFTruth *myDHDDMTOFTruth = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDHDDMTOFTruth->x);
	//			printcol("%3.2f",	myDHDDMTOFTruth->y);
	//			printrow();
	//		}
	//
	return _table;

}




//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DHDDMTOFTruth::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{

	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
    v.clear();
	
        // Loop over Physics Events
    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE) return NOERROR;
	
    identifier_t identifier = 0;
	
    for(unsigned int i=0; i<PE->mult; i++){

        s_TofPoints_t*  tofpoints = NULL;
        if(PE->in[i].hitView)
            if(PE->in[i].hitView->forwardTOF)
                tofpoints = PE->in[i].hitView->forwardTOF->tofPoints;
        if(!tofpoints)continue;
        
        for(unsigned int j=0;j<tofpoints->mult;j++){

            DHDDMTOFTruth *toftruth = new DHDDMTOFTruth;
            toftruth->orientation = 0;
            toftruth->track       = tofpoints->in[j].track;
            toftruth->primary     = tofpoints->in[j].primary ? 1 : 0;
            toftruth->x           = tofpoints->in[j].x;
            toftruth->y           = tofpoints->in[j].y;
            toftruth->z           = tofpoints->in[j].z;
            toftruth->t           = tofpoints->in[j].t;
            identifier++;
            v.push_back(toftruth);

        }
    }

    return NOERROR;

}

