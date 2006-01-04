// $Id$
//
//    File: DFactory_DFCALTruthShower.cc
// Created: Wed Jan  4 14:43:05 EST 2006
// Creator: davidl (on Linux jlabl1.jlab.org 2.4.21-37.ELsmp i686)
//

#include <cassert>	

#include "DFactory_DFCALTruthShower.h"

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DFCALTruthShower::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardEMcal == HDDM_NULL ||
			hits->forwardEMcal->fcalTruthShowers == HDDM_NULL)continue;

		s_FcalTruthShowers_t *showers = hits->forwardEMcal->fcalTruthShowers;
		for(unsigned int j=0; j<showers->mult; j++){
			s_FcalTruthShower_t *shower = &showers->in[j];
			
			DFCALTruthShower *dfcaltruthshower = new DFCALTruthShower(
				shower->x,
				shower->y,
				shower->z,
				shower->E,
				shower->t,
				shower->primary,
				shower->track
				);
			
			_data.push_back(dfcaltruthshower);
		}

	} // i  (physicsEvents)

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DFCALTruthShower::toString(void)
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
	//			DFCALTruthShower *myDFCALTruthShower = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDFCALTruthShower->x);
	//			printcol("%3.2f",	myDFCALTruthShower->y);
	//			printrow();
	//		}
	//
	return _table;

}
