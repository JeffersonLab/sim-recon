// $Id$
//
//    File: DFactory_DTOFPoint.cc
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>	

#include "DFactory_DTOFPoint.h"
#include "DFactory_DTOFHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DTOFPoint::evnt(DEventLoop *loop, int eventnumber)
{

  assert( _data.size() == 0 );

  vector<const DTOFHit*> hits;
  eventLoop->Get(hits);

  for (unsigned int i = 0; i < hits.size(); i++){

    const DTOFHit *hit = hits[i];
    DTOFPoint *point = new DTOFPoint;

    point->x         = 0.0;
    point->y         = hit->y;
    point->z         = 0.0;
    point->E         = hit->E;
    point->t         = hit->t;
    point->nhits     = 1;
    point->hits[0]   = i;
    point->chisq     = 0.0;

    _data.push_back(point);

  }

  return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DTOFPoint::toString(void)
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
	//			DTOFPoint *myDTOFPoint = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDTOFPoint->x);
	//			printcol("%3.2f",	myDTOFPoint->y);
	//			printrow();
	//		}
	//
	return _table;

}
