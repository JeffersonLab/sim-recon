// $Id$
//
//    File: DTOFPoint_factory.cc
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>
#include <math.h>
using namespace std;

#include "DTOFPoint_factory.h"
#include "DTOFTruth.h"
#include "DTOFHit.h"

#define MAXTOFHITS 50
#define VELOCITY 15.0
#define LONGBARLENGTH 258.0
#define BARWIDTH 6.0

//------------------
// evnt
//------------------
jerror_t DTOFPoint_factory::evnt(JEventLoop *loop, int eventnumber)
{

  assert( _data.size() == 0 );

  vector<const DTOFHit*> hits;
  eventLoop->Get(hits);

  vector<const DTOFTruth*> tracks;
  eventLoop->Get(tracks);

  unsigned int usedlist[MAXTOFHITS*4];
  for (unsigned int i = 0; i < hits.size(); i++){
  	if(i>=MAXTOFHITS*4){cout<<__FILE__<<":"<<__LINE__<<" too many hits in TOF!!"<<endl; break;}
    usedlist[i] = 0;
  }

  for (unsigned int i = 0; i < tracks.size(); i++){
    float x = tracks[i]->x;
    float y = tracks[i]->y;
    float z = tracks[i]->z;
    unsigned int trackid = i;
    unsigned int nhits = 0;
    unsigned int hitlist[16];
    int ptype = tracks[i]->ptype;
    for (unsigned int j = 0; j < hits.size(); j++){
      //      float y0 = hits[j]->y;
      // y does not exist anymore this bas to be handled differently in the future
      // by using the bar number hits[j]->bar
      float y0 = (hits[j]->bar-20.5)*6.; // here the half bars are not taken into account!!!!!

      if (((hits[j]->orientation == 0)&&(fabs(x-y0)<=1.0*BARWIDTH))||
          ((hits[j]->orientation == 1)&&(fabs(y-y0)<=1.0*BARWIDTH))){
		  if(nhits>=16){cout<<__FILE__<<":"<<__LINE__<<" too many hits in TOF!!"<<endl; break;}
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
      point->ptype = ptype;
      _data.push_back(point);
    }
  }

  return NOERROR;
}

