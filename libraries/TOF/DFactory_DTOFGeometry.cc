// $Id$
//
//    File: DFactory_DTOFGeometry.cc
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#include "DFactory_DTOFGeometry.h"

//------------------
// evnt
//------------------
derror_t DFactory_DTOFGeometry::evnt(DEventLoop *loop, int eventnumber)
{
	// Code to generate factory data goes here. Add it like:
	//
	// DTOFGeometry *myDTOFGeometry = new DTOFGeometry;
	// myDTOFGeometry->x = x;
	// myDTOFGeometry->y = y;
	// ...
	// _data.push_back(myDTOFGeometry);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DTOFGeometry::toString(void)
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
	//			DTOFGeometry *myDTOFGeometry = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDTOFGeometry->x);
	//			printcol("%3.2f",	myDTOFGeometry->y);
	//			printrow();
	//		}
	//
	return _table;

}
