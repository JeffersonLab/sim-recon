// $Id$
//
//    File: DFCALGeometry_factory.cc
// Created: Wed Aug 24 10:09:27 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#include <cassert>

#include "DFCALGeometry_factory.h"
#include "DFCALGeometry.h"

//------------------
// brun
//------------------
jerror_t DFCALGeometry_factory::brun(JEventLoop *loop, int runnumber)
{
	assert( _data.size() == 0 );

	flags = PERSISTANT;
	_data.push_back( new DFCALGeometry() );
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DFCALGeometry_factory::erun(JEventLoop *loop)
{
	for(unsigned int i=0; i<_data.size(); i++)delete _data[i];
	_data.clear();
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFCALGeometry_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The JFactory_base class has several methods defined to help. They
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
