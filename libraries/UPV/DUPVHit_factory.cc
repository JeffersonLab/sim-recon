// $Id$
//
//    File: DUPVHit_factory.cc
// Created: Thu Jun  9 10:01:38 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DUPVHit_factory.h"

//------------------
// Etotal
//------------------
double DUPVHit_factory::Etotal(void)
{
	// Assume the objects have already been retrieved
	double E=0.0;
	for(unsigned int i=0; _data.size(); i++)E+=_data[i]->E;
	
	return E;
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
	printheader("layer:    row:      E(MeV):      t(ns):    side:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DUPVHit *hit = _data[i];
	
		printnewrow();
		printcol("%d",	hit->layer);
		printcol("%d",	hit->row);
		printcol("%1.4f",	hit->E*1000.0);
		printcol("%3.2f",	hit->t);
		printcol("%s",	hit->side==DUPVHit::UPV_LEFT ? "left":"right");
		printrow();
	}

	return _table;
}
