// $Id$
//
//    File: DBCALHit_factory.cc
// Created: Thu Jun  9 10:14:35 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DBCALHit_factory.h"

//------------------
// toString
//------------------
const string DBCALHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()==0)return string(); // don't print anything if we have no data!

	printheader("row:   module:  layer:  sector:         end:     E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DBCALHit *bcalhit = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%d",	bcalhit->module);
		printcol("%d",	bcalhit->layer);
		printcol("%d",	bcalhit->sector);
		printcol(bcalhit->end==DBCALHit::UPSTREAM ? "upstream":"downstream");
		printcol("%2.3f",	bcalhit->E);
		printcol("%4.0f",	bcalhit->t);
		printrow();
	}

	return _table;
}
