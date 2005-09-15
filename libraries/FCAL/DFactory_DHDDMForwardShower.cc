// $Id$
//
//    File: DFactory_DHDDMForwardShower.cc
// Created: Mon Aug 29 15:14:08 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#include "DFactory_DHDDMForwardShower.h"

//------------------
// evnt
//------------------
derror_t DFactory_DHDDMForwardShower::evnt(DEventLoop *loop, int eventnumber)
{
	// no code should be here -- this factory is used strictly for reading in
	// HDDM data
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DHDDMForwardShower::toString(void)
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
	//			DHDDMForwardShower *myDHDDMForwardShower = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDHDDMForwardShower->x);
	//			printcol("%3.2f",	myDHDDMForwardShower->y);
	//			printrow();
	//		}
	//
	return _table;

}

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DHDDMForwardShower::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
	v.clear();
	
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	identifier_t identifier = 0;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_Rows_t *rows = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardEMcal)
				rows = PE->in[i].hitView->forwardEMcal->rows;
		if(!rows)continue;
		
		for(unsigned int j=0;j<rows->mult;j++){
			float y = rows->in[j].y;
			s_Columns_t *columns = rows->in[j].columns;
			if(!columns)continue;
			
			for(unsigned int k=0;k<columns->mult;k++){
				float x = columns->in[k].x;
				s_Showers_t *showers = columns->in[k].showers;
				if(!showers)continue;
				
				for(unsigned int m=0;m<showers->mult;m++){
					float E = showers->in[m].E;
					float t = showers->in[m].t;
					
					v.push_back( new DHDDMForwardShower( identifier++, 
														 x, y, E, t ) );
				}
			}
		}
	}
	
	return NOERROR;
}
