// $Id$
//
//    File: DUPVTruthHit_factory.cc
// Created: Mon Nov  6 09:58:35 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#include "DUPVTruthHit_factory.h"

//------------------
// GetETotal
//------------------
double DUPVTruthHit_factory::GetETotal(void)
{
	Get();
	
	double Etotal = 0.0;
	for(unsigned int i=0; i<_data.size(); i++)Etotal+=_data[i]->E;
	
	return Etotal;
}


//------------------
// toString
//------------------
const string DUPVTruthHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("  E(MeV): primary:   t (ns): track:       x:       y:       z:");

	for(unsigned int i=0; i<_data.size(); i++){
		DUPVTruthHit *hit = _data[i];

		printnewrow();
		printcol("%3.3f",	hit->E*1000.0);
		printcol("%d",	hit->primary);
		printcol("%2.2f",	hit->t);
		printcol("%d",	hit->track);
		printcol("%3.1f",	hit->x);
		printcol("%3.1f",	hit->y);
		printcol("%3.1f",	hit->z);
		printrow();
	}

	return _table;
}
