// $Id$
//
//    File: DFactory_DUPVHit.cc
// Created: Sun Apr  3 10:35:31 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DFactory_DUPVHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DUPVHit::evnt(int enventnumber)
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
const string DFactory_DUPVHit::toString(void)
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
	//			DUPVHit *myDUPVHit = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDUPVHit->x);
	//			printcol("%3.2f",	myDUPVHit->y);
	//			printrow();
	//		}
	//
	//		return _table;

}
