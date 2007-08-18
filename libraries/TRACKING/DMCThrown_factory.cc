// $Id$
//
//    File: DMCThrown_factory.cc
// Created: Sun Apr  3 12:22:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DMCThrown_factory.h"

//------------------
// evnt
//------------------
jerror_t DMCThrown_factory::evnt(JEventLoop *loop, int enventnumber)
{
	/// This doesn't do anything. All of the objects are obtained
	/// from the event source. This class exists only to to provide 
	/// a toString method for pretty printing.

	return NOERROR;
}

//------------------
// toString
//------------------
const string DMCThrown_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("  myid: parent: type: pdgtype: mech:  q:     p:    E: theta:   phi:   mass:     x:     y:     z:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DMCThrown * mcthrown = _data[i];

		printnewrow();
		
		printcol("%d", mcthrown->myid);
		printcol("%d", mcthrown->parentid);
		printcol("%d", mcthrown->type);
		printcol("%d", mcthrown->pdgtype);
		printcol("%d", mcthrown->mech);
		printcol("%+d", (int)mcthrown->q);
		printcol("%3.3f", mcthrown->p);
		printcol("%3.1f", mcthrown->E);
		printcol("%1.3f", mcthrown->theta);
		printcol("%1.3f", mcthrown->phi);
		printcol("%1.3f", mcthrown->mass);
		printcol("%2.2f", mcthrown->x);
		printcol("%2.2f", mcthrown->y);
		printcol("%2.2f", mcthrown->z);

		printrow();
	}
	
	return _table;
}
