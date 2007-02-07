// $Id$
//
//    File: DFCALHit_factory.cc
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <cassert>

#include "DFCALHit_factory.h"
#include "DFCALHit.h"
#include "DFCALMCResponse.h"

//------------------
// toString
//------------------
const string DFCALHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("hit:  column:   row:   x(cm):   y(cm):   E(MeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALHit *fcalhit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%d", fcalhit->column);
		printcol("%d", fcalhit->row);
		printcol("%3.1f", fcalhit->x);
		printcol("%3.1f", fcalhit->y);
		printcol("%2.3f", fcalhit->E*1000.0);
		printcol("%4.0f", fcalhit->t);
		printrow();
	}
	
	return _table;
}
