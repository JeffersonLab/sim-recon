// $Id$
//
//    File: DSCTruthHit_factory.cc
// Created: Wed Feb  7 10:53:46 EST 2007
// Creator: davidl (on Linux megrez.jlab.org 2.6.9-42.0.2.ELsmp x86_64)
//


#include "DSCTruthHit_factory.h"


//------------------
// toString
//------------------
const string DSCTruthHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ADSCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	printheader("track: primary:    dEdx(MeV/cm):      t:      r:     phi:     z:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DSCTruthHit *hit = _data[i];
	
		printnewrow();
		printcol("%d",	hit->track);
		printcol("%d",	hit->primary);
		printcol("%1.3f",	hit->dEdx*1.0E3);
		printcol("%3.2f",	hit->t);
		printcol("%3.1f",	hit->r);
		printcol("%1.3f",	hit->phi);
		printcol("%3.1f",	hit->z);
		printrow();
	}

	return _table;
}
