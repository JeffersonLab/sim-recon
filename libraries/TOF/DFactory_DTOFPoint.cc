// $Id$
//
//    File: DFactory_DTOFPoint.cc
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>	

#include "DFactory_DHDDMTOFTruth.h"
#include "DFactory_DTOFPoint.h"
#include "DFactory_DTOFHit.h"

#define MAXTOFHITS 50
#define VELOCITY 15.0
#define LONGBARLENGTH 258.0
#define BARWIDTH 6.0

//------------------
// evnt
//------------------
derror_t DFactory_DTOFPoint::evnt(DEventLoop *loop, int eventnumber)
{

  assert( _data.size() == 0 );

  vector<const DTOFHit*> hits;
  eventLoop->Get(hits);

  vector<const DHDDMTOFTruth*> tracks;
  eventLoop->Get(tracks);

  unsigned int usedlist[MAXTOFHITS*4];
  for (unsigned int i = 0; i < hits.size(); i++){
    usedlist[i] = 0;
  }

  for (unsigned int i = 0; i < tracks.size(); i++){
    float x = tracks[i]->x;
    float y = tracks[i]->y;
    float z = tracks[i]->z;
    unsigned int trackid = i;
    unsigned int nhits = 0;
    unsigned int hitlist[16];
    for (unsigned int j = 0; j < hits.size(); j++){
      float y0 = hits[j]->y;
      if (((hits[j]->orientation == 0)&&(fabs(x-y0)<=1.0*BARWIDTH))||
          ((hits[j]->orientation == 1)&&(fabs(y-y0)<=1.0*BARWIDTH))){
        hitlist[nhits] = j;
        nhits++;
      }
    }
    if (nhits > 0){
      DTOFPoint *point = new DTOFPoint;
      point->x = x;
      point->y = y;
      point->z = z;
      point->nhits = nhits;
      for (unsigned int j = 0; j < nhits; j++){point->hits[j] = hitlist[j];}
      point->dedx = 0.0;
      point->chisq = 0.0;
      point->t = 0.0;
      point->trackid = trackid;
      _data.push_back(point);
    }
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
