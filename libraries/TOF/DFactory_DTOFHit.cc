// $Id$
//
//    File: DFactory_DTOFHit.cc
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DFactory_DTOFHit.h"
#include "DFactory_DTOFMCResponse.h"

//------------------
// evnt
//------------------
derror_t DFactory_DTOFHit::evnt(DEventLoop *eventLoop, int eventnumber)
{

  vector<const DTOFMCResponse*> mcresponses;
  eventLoop->Get(mcresponses);

  for (unsigned int i = 0; i < mcresponses.size(); i++){

    const DTOFMCResponse *mcresponse = mcresponses[i];
    DTOFHit *hit = new DTOFHit;

    hit->orientation = mcresponse->orientation;
    hit->end         = mcresponse->end;
    hit->y           = mcresponse->y;
    hit->t           = mcresponse->t;
    hit->E           = mcresponse->E;

    _data.push_back(hit);

  }

  return NOERROR;
}


//------------------
// toString
//------------------
const string DFactory_DTOFHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   x(cm):   y(cm):  orientation:     end:     dE(MeV):   t(ns):");

/*	
	for(unsigned int i=0; i<_data.size(); i++){
		DTOFHit *tofhit = _data[i];

		printnewrow();
		
		printcol("%d", i);
		printcol("%3.1f", tofhit->x);
		printcol("%3.1f", tofhit->y);
		printcol(tofhit->orientation ? "horizontal":"vertical");
		printcol(tofhit->end ? "right":"left");
		printcol("%2.3f", tofhit->dE*1000.0);
		printcol("%4.3f", tofhit->t);
		
		printrow();
	}
*/
	
	return _table;
}
