// $Id$
//
//    File: DFactory_DTOFMCResponse.cc
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#include "DFactory_DTOFMCResponse.h"
#include "DFactory_DHDDMTOFHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DTOFMCResponse::evnt(DEventLoop *loop, int eventnumber)
{

  vector<const DHDDMTOFHit*> hddmhits;
  eventLoop->Get(hddmhits);

  for (unsigned int i = 0; i < hddmhits.size(); i++){

    const DHDDMTOFHit *hddmhit = hddmhits[i];
    DTOFMCResponse *response = new DTOFMCResponse;

    response->orientation = hddmhit->orientation;
    response->end         = hddmhit->end;
    response->y           = hddmhit->y;
    response->t           = hddmhit->t;
    response->E           = hddmhit->E;

    _data.push_back(response);

  }

  return NOERROR;
}


//------------------
// toString
//------------------
const string DFactory_DTOFMCResponse::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	//		printheader("row:    x:     y:");
	//
	// 	for(int i=0; i<_data.size(); i++){
	//			DTOFMCResponse *myDTOFMCResponse = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDTOFMCResponse->x);
	//			printcol("%3.2f",	myDTOFMCResponse->y);
	//			printrow();
	//		}
	//
	return _table;

}
