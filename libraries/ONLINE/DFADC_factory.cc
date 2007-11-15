// $Id$
//
//    File: DFADC_factory.cc
// Created: Thu Nov 15 09:51:29 EST 2007
// Creator: bellis (on Linux mordor 2.6.22.1 unknown)
//


#include "DFADC_factory.h"

//------------------
// evnt
//------------------
jerror_t DFADC_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// DFADC *myDFADC = new DFADC;
	// myDFADC->x = x;
	// myDFADC->y = y;
	// ...
	// _data.push_back(myDFADC);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFADC_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	printheader("row:    x:     y:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFADC *myDFADC = _data[i];
	
		printnewrow();
		printcol("%d",	i);
//		printcol("%1.3f",	myDFADC->x);
//		printcol("%3.2f",	myDFADC->y);
		printrow();
	}

	return _table;
}
