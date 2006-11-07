// $Id$
//
//    File: DFCALTruthShower_factory.cc
// Created: Wed Jan  4 14:43:05 EST 2006
// Creator: davidl (on Linux jlabl1.jlab.org 2.4.21-37.ELsmp i686)
//

#include <cassert>	

#include "DFCALTruthShower_factory.h"

//------------------
// toString
//------------------
const string DFCALTruthShower_factory::toString(void)
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
