// $Id$
//
//    File: DCDCHit_factory.cc
// Created: Thu Jun  9 10:22:37 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DCDCHit_factory.h"


//------------------
// toString
//------------------
const string DCDCHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: ring:  straw:  dE(MeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DCDCHit *cdchit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%d", cdchit->ring);
		printcol("%d", cdchit->straw);
		printcol("%2.3f", cdchit->dE);
		printcol("%4.0f", cdchit->t);
		printrow();
	}

	return _table;
}
