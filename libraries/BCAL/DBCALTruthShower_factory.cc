// $Id$
//
//    File: DBCALTruthShower_factory.cc
// Created: Fri Nov 18 10:37:37 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#include <cassert>	

#include "DBCALTruthShower_factory.h"


//------------------
// toString
//------------------
const string DBCALTruthShower_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("track:   phi:          r:          z:        E:        t:");

	for(unsigned int i = 0; i < _data.size(); i++) {
		DBCALTruthShower *s = _data[i];

		printnewrow();
		printcol("%i",s->track);
		printcol("%4.3f",s->phi);
		printcol("%3.1f",s->r);
		printcol("%4.1f",s->z);
		printcol("%4.3f",s->E);
		printcol("%4.3f",s->t);
		printrow();
	}

  printnewrow();
  printrow();

     return _table;

}


