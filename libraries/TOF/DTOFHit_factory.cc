// $Id$
//
//    File: DTOFHit_factory.cc
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DTOFHit_factory.h"
#include "DTOFMCResponse.h"

//------------------
// evnt
//------------------
jerror_t DTOFHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{

  vector<const DTOFMCResponse*> mcresponses;
  eventLoop->Get(mcresponses);

  for (unsigned int i = 0; i < mcresponses.size(); i++){

    const DTOFMCResponse *mcresponse = mcresponses[i];
    DTOFHit *hit = new DTOFHit;

    // do any run dependent calibration here

    hit->id          = mcresponse->id;
    hit->orientation = mcresponse->orientation;
    hit->t_north     = mcresponse->t_north;
    hit->E_north     = mcresponse->E_north;
    hit->t_south     = mcresponse->t_south;
    hit->E_south     = mcresponse->E_south;

    _data.push_back(hit);

  }

  return NOERROR;
}


//------------------
// toString
//------------------
const string DTOFHit_factory::toString(void)
{
  // Ensure our Get method has been called so _data is up to date
  Get();
  if(_data.size()<=0)return string(); // don't print anything if we have no data!

  printheader( "id: orientation: end:    t [ns]:    x/y (orth.):   dE [MeV]:" );

	
  for(unsigned int i=0; i<_data.size(); i++){
    DTOFHit *tofhit = _data[i];

    printnewrow();
    printcol("%d",	tofhit->id );
    printcol("%d",	tofhit->orientation );
    printcol("%1.3f",	tofhit->t_north );
    printcol("%1.3f",	tofhit->E_north );    
    printcol("%1.3f",	tofhit->t_south );
    printcol("%1.3f",	tofhit->E_south );    

    printrow();
  }
  
	
  return _table;
}
