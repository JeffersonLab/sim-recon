// $Id$
//
//    File: DSCHit_factory.cc
// Created: Wed Feb  7 10:46:20 EST 2007
// Creator: davidl (on Linux megrez.jlab.org 2.6.9-42.0.2.ELsmp x86_64)
//


#include "DSCHit_factory.h"


//------------------
// toString
//------------------
const string DSCHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ADSCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	printheader("sector:    t:   dE(MeV):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DSCHit *hit = _data[i];
	
		printnewrow();
		printcol("%d",	hit->sector);
		printcol("%2.3f",	hit->t);
		printcol("%1.3f",	hit->dE*1.0E3);
		printrow();
	}

	return _table;
}
