// $Id$
//
//    File: DTOFMCResponse_factory.cc
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#include "DTOFMCResponse_factory.h"
#include "DHDDMTOFHit.h"
#include "DTOFGeometry.h"

//------------------
// evnt
//------------------
jerror_t DTOFMCResponse_factory::evnt(JEventLoop *loop, int eventnumber)
{

  vector<const DHDDMTOFHit*> hddmhits;
  eventLoop->Get(hddmhits);

  vector<const DTOFGeometry*> tofGeomVect;
  eventLoop->Get(tofGeomVect);
 
  const DTOFGeometry& tofGeom = (*(tofGeomVect[0]));

  for (unsigned int i = 0; i < hddmhits.size(); i++){

    const DHDDMTOFHit *hddmhit = hddmhits[i];
    DTOFMCResponse *response = new DTOFMCResponse;

    // do any run-dependent smearing here

    response->id          = hddmhit->id;
    response->orientation = hddmhit->plane;
    response->end         = hddmhit->end;
    response->y           = tofGeom.bar2y( hddmhit->bar, hddmhit->plane );
    response->t           = hddmhit->t;
    response->E           = hddmhit->dE;

    _data.push_back(response);

  }

  return NOERROR;
}


//------------------
// toString
//------------------
const string DTOFMCResponse_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
  Get();
  if(_data.size()<=0)return string(); // don't print anything if we have no data!

  // Put the class specific code to produce nicely formatted ASCII here.
  // The JFactory_base class has several methods defined to help. They
  // rely on positions of colons (:) in the header. Here's an example:
  
  printheader( "id: orientation: end:    t [ns]:    x/y (orth.):   dE [MeV]:" );

  for(unsigned int i=0; i<_data.size(); i++ ){

    DTOFMCResponse *myTOF = _data[i];
    
    printnewrow();
    printcol("%d",	myTOF->id );
    printcol("%d",	myTOF->orientation );
    printcol("%d",	myTOF->end );
    printcol("%1.3f",	myTOF->t );
    printcol("%2.3f",	myTOF->y );
    printcol("%1.3f",	myTOF->E );
    printrow();
  }
  
  return _table;

}
