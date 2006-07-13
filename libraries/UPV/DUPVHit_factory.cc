// $Id$
//
//    File: DUPVHit_factory.cc
// Created: Thu Jun  9 10:01:38 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DUPVHit_factory.h"

//------------------
// evnt
//------------------
jerror_t DUPVHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	// Code to generate factory data goes here. Add it like:
	//
	// DUPVHit *myDUPVHit = new DUPVHit;
	// myDUPVHit->x = x;
	// myDUPVHit->y = y;
	// ...
	// _data.push_back(myDUPVHit);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// toString
//------------------
const string DUPVHit_factory::toString(void)
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
	//			DUPVHit *myDUPVHit = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDUPVHit->x);
	//			printcol("%3.2f",	myDUPVHit->y);
	//			printrow();
	//		}
	//
	return _table;

}
