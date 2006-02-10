// $Id$
//
//    File: DFactory_DHDDMBCALTruth.cc
// Created: Fri Nov 18 10:37:37 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#include <cassert>	

#include "DFactory_DHDDMBCALTruth.h"

//------------------
// evnt
//------------------
derror_t DFactory_DHDDMBCALTruth::evnt(DEventLoop *loop, int eventnumber)
{
  // no code should be here -- this factory is used strictly for reading in
  // HDDM data

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DHDDMBCALTruth::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

  printheader("#:   phi:          r:          z:        E:        t:");

  for(unsigned int i = 0; i < _data.size(); i++) {
    DHDDMBCALTruth *s = _data[i];

/////////////////////////////////////////////
    printf("%d  %f  %f  %f  %f  %f\n",s->track,s->phi,s->r,s->z,s->E,s->t);
// Commented out until all core dumps in this printing method are gone
//    printnewrow();
//    printcol("%i",s->track);
//    printcol("%f",s->phi);
//    printcol("%f",s->r);
//    printcol("%f",s->z);
//    printcol("%f",s->E);
//    printcol("%f",s->t);
//    printrow();
  }

  printnewrow();
  printrow();

     return _table;

}


//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DHDDMBCALTruth::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
  // Copies the data from the given hddm_s structure. This is called
  // from DEventSourceHDDM::GetObjects.

  v.clear();
	
  // Loop over Physics Events
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) 
    return NOERROR;
	
  identifier_t identifier = 0;
	
  for(unsigned int i = 0; i < PE->mult; i++) {
    s_BcalTruthShowers_t* bcalTruthShowers = NULL;
    if(PE->in[i].hitView)
      if(PE->in[i].hitView->barrelEMcal)
        bcalTruthShowers= PE->in[i].hitView->barrelEMcal->bcalTruthShowers;
    if(!bcalTruthShowers)continue;
        
    for(unsigned int j = 0; j < bcalTruthShowers->mult; j++) {
      DHDDMBCALTruth *bcaltruth = new DHDDMBCALTruth;
      bcaltruth->track = bcalTruthShowers->in[j].track;
      bcaltruth->primary = bcalTruthShowers->in[j].primary ? 1 : 0;
      bcaltruth->phi = bcalTruthShowers->in[j].phi;
      bcaltruth->r = bcalTruthShowers->in[j].r;
      bcaltruth->z = bcalTruthShowers->in[j].z;
      bcaltruth->t = bcalTruthShowers->in[j].t;
      bcaltruth->E = bcalTruthShowers->in[j].E;
      identifier++;
      v.push_back(bcaltruth);
    }
  }

  return NOERROR;
}
