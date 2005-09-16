// $Id$
//
//    File: DFactory_DFCALGeometry.cc
// Created: Wed Aug 24 10:09:27 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#include <cassert>

#include "DFactory_DFCALGeometry.h"
#include "DFCALGeometry.h"

//------------------
// evnt
//------------------
derror_t DFactory_DFCALGeometry::evnt(DEventLoop *loop, int eventnumber)
{

	assert( _data.size() == 0 );
		
	_data.push_back( new DFCALGeometry() );
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DFCALGeometry::toString(void)
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
	//			DFCALGeometry *myDFCALGeometry = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDFCALGeometry->x);
	//			printcol("%3.2f",	myDFCALGeometry->y);
	//			printrow();
	//		}
	//
	return _table;

}
