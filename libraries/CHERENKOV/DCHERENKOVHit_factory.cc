// $Id$
//
//    File: DCHERENKOVHit_factory.cc
// Created: Thu Jun  9 10:32:49 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DCHERENKOVHit_factory.h"

//------------------
// evnt
//------------------
jerror_t DCHERENKOVHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	/// Place holder for now. 

	return NOERROR;
}

//------------------
// Extract_HDDM
//------------------
jerror_t DCHERENKOVHit_factory::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// THIS NEEDS TO BE WRITTEN!!
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects.
	
	v.clear();

	return NOERROR;
}

//------------------
// toString
//------------------
const string DCHERENKOVHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	// GetNrows() will check the data source first in case the objects
	// are obtained fom there.
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The JFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	//		printheader("row:    x:     y:");
	//
	// 	for(int i=0; i<_data.size(); i++){
	//			DCHERENKOVHit *myDCHERENKOVHit = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDCHERENKOVHit->x);
	//			printcol("%3.2f",	myDCHERENKOVHit->y);
	//			printrow();
	//		}
	//
	return _table;

}
