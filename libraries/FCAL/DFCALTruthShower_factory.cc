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

	printheader("id:      x:      y:      z:   E(MeV): t(ns): primary: track:");

	for(unsigned int i=0; i<_data.size(); i++){
		DFCALTruthShower *shower = _data[i];

		printnewrow();
		printcol("%d",	shower->id);
		printcol("%3.1f",	shower->x());
		printcol("%3.1f",	shower->y());
		printcol("%3.1f",	shower->z());
		printcol("%3.3f",	shower->E()*1000.0);
		printcol("%3.1f",	shower->t());
		printcol("%d",	shower->primary());
		printcol("%d",	shower->track());
		printrow();
	}

	return _table;
}
