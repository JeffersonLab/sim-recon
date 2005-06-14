// $Id$
//
//    File: DFactory_DTRIGGER.cc
// Created: Thu Jun  9 10:03:17 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DFactory_DTRIGGER.h"

//------------------
// evnt
//------------------
derror_t DFactory_DTRIGGER::evnt(DEventLoop *eventLoop, int eventnumber)
{
	// Code to generate factory data goes here. Add it like:
	//
	// DTRIGGER *myDTRIGGER = new DTRIGGER;
	// myDTRIGGER->x = x;
	// myDTRIGGER->y = y;
	// ...
	// _data.push_back(myDTRIGGER);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DTRIGGER::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	// GetNrows() will check the data source first in case the objects
	// are obtained fom there.
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	//		printheader("row:    x:     y:");
	//
	// 	for(int i=0; i<_data.size(); i++){
	//			DTRIGGER *myDTRIGGER = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDTRIGGER->x);
	//			printcol("%3.2f",	myDTRIGGER->y);
	//			printrow();
	//		}
	//
	return _table;

}
