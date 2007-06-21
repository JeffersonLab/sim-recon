// $Id$
//
//    File: DHDDMTOFHit_factory.cc
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//
// changes: Tue Jun 19 17:29:07 EDT 2007 B.Zihlmann
//          adjust that the hit structure contains both north and south 
//          information as well as the hit location.
//

#include <cassert>	

#include "DHDDMTOFHit_factory.h"

//------------------
// evnt
//------------------
jerror_t DHDDMTOFHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// no code should be here -- this factory is used strictly for reading in
	// HDDM data
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DHDDMTOFHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("id:  bar:  plane:   end:        t:     dE:   ");

	for(unsigned int i=0; i<_data.size(); i++){
       	  DHDDMTOFHit *hit = _data[i];

		printnewrow();
		printcol("%d",		hit->id);
		printcol("%d",		hit->bar);
		printcol("%d",		hit->plane);
		printcol("%1.3f",	hit->t_north);
		printcol("%1.3f",	hit->dE_north);
		printcol("%1.3f",	hit->t_south);
		printcol("%1.3f",	hit->dE_south);
		printcol("%12.4e",	hit->x);
		printcol("%12.4e",	hit->y);
		printcol("%12.4e",	hit->z);
		printcol("%12.4e",	hit->px);
		printcol("%12.4e",	hit->py);
		printcol("%12.4e",	hit->pz);
		printcol("%12.4e",	hit->E);
		printcol("%12.4e",	hit->ptype);
		printrow();
	}

	return _table;
}


