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

#define MAXTOFHITS 50
#define VELOCITY 15.0
#define LONGBARLENGTH 252.0
#define BARWIDTH 6.0


//------------------
// evnt
//------------------
jerror_t DTOFPoint_factory::evnt(JEventLoop *loop, int eventnumber)
{
  // Get TOF hits
  vector<const DHDDMTOFHit*> hits;
  eventLoop->Get(hits);

  // Local vectors for separating hits by plane
  vector<const DHDDMTOFHit *>uhits;
  vector<const DHDDMTOFHit *>vhits;
  for (unsigned int j = 0; j < hits.size(); j++){
    const DHDDMTOFHit *hit=hits[j];

    if (hit->plane) uhits.push_back(hit);
    else vhits.push_back(hit);
  }

  // Match hits in the two planes by bar number and position as determined from
  // the time difference between two ends, where available.
  for (unsigned int i=0;i<uhits.size();i++){
    double ux=VELOCITY*(uhits[i]->t_south-uhits[i]->t_north)/2.;
    double utof=(uhits[i]->t_north+uhits[i]->t_south)/2.
      -LONGBARLENGTH/2./VELOCITY;
    double dy=uhits[i]->bar>20?BARWIDTH:-BARWIDTH;
    double uy=BARWIDTH*(uhits[i]->bar-20.5)+dy;
    double x_cut=BARWIDTH*0.55; // changed to "*0.55" from "/2." 4/28/2010  DL
    double t_cut=1.;
    
    // Handle the single-ended bars
    if (uhits[i]->bar>40){
      x_cut=120.;
      t_cut=10.;
      switch(uhits[i]->bar){
      case 41:
	uy=-3.;
	ux=-72.;
	utof=uhits[i]->t_south-LONGBARLENGTH/2./VELOCITY;
	break;
      case 43:
	uy=3.;
	ux=-72.;
	utof=uhits[i]->t_south-LONGBARLENGTH/2./VELOCITY;
	break;
      case 42:
	uy=-3.;	
	ux=+72.;
	utof=uhits[i]->t_north-LONGBARLENGTH/2./VELOCITY;
	break;
      case 44:
	uy=+3.;
	ux=+72.;
	utof=uhits[i]->t_north-LONGBARLENGTH/2./VELOCITY;
	break;
      default:
	break;
      }

    }

    for (unsigned int j=0;j<vhits.size();j++){
      double vy=VELOCITY*(vhits[j]->t_south-vhits[j]->t_north)/2.;
      double vtof=(vhits[j]->t_north+vhits[j]->t_south)/2.
	-LONGBARLENGTH/2./VELOCITY;
      double dx=vhits[j]->bar>20?BARWIDTH:-BARWIDTH;
      double vx=BARWIDTH*(vhits[j]->bar-20.5)+dx;
      double y_cut=BARWIDTH*0.55; // changed to "*0.55" from "/2." 4/28/2010  DL

      // Handle the single-ended bars
      if (vhits[j]->bar>40){
	y_cut=120.;
	t_cut=10.;
	switch(vhits[j]->bar){
	case 41:
	  vx=-3.;
	  vy=72.;
	  vtof=vhits[j]->t_south-LONGBARLENGTH/2./VELOCITY;
	  break;
	case 43:
	  vx=3.;
	  vy=72.;
	  vtof=vhits[j]->t_south-LONGBARLENGTH/2./VELOCITY;
	  break;
	case 42:
	  vx=-3.;
	  vy=-72.;
	  vtof=vhits[j]->t_north-LONGBARLENGTH/2./VELOCITY;
	  break;
	case 44:
	  vx=3.;
	  vy=-72;
	  vtof=vhits[j]->t_north-LONGBARLENGTH/2./VELOCITY;
	  break;
	default:
	  break;
	}	
      }
      
      // If we have a match in both position and time, output the point with 
      // the combined 2-layer data
      if (fabs(ux-vx)<x_cut && fabs(uy-vy)<y_cut && fabs(utof-vtof)<t_cut){
	 DTOFPoint *point = new DTOFPoint;

	 point->chisq=0.; //??
	 if (vhits[j]->bar<41 && uhits[i]->bar<41){
	   point->pos(1)=vy;
	   point->pos(0)=ux;
	   point->pos(2)=618.81;
	   point->t=(utof+vtof)/2.;
	   point->dedx=(sqrt(uhits[i]->dE_north*uhits[i]->dE_south)
			+sqrt(vhits[j]->dE_north*vhits[j]->dE_south))/2.;
	 }
	 else if (vhits[j]->bar<41 && uhits[i]->bar>40){   
	   point->pos(0)=vx;
	   point->pos(1)=vy;
	   point->pos(2)=617.52;
	   point->t=vtof;
	   point->dedx=sqrt(vhits[j]->dE_north*vhits[j]->dE_south);
	 }
	 else{
	   point->pos(0)=ux;
	   point->pos(1)=uy;
	   point->pos(2)=620.10;
	   point->t=utof;
	   point->dedx=sqrt(uhits[i]->dE_north*uhits[i]->dE_south);
	 }
      	 _data.push_back(point);
      }
    }
  }
  
  return NOERROR;
}

