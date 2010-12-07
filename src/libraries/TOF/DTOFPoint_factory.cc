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

//------------------
// brun
//------------------
jerror_t DTOFPoint_factory::brun(JEventLoop *loop, int runnumber)
{

  map<string, double> tofparms;
 
  if ( !loop->GetCalib("TOF/tof_parms", tofparms)){
    cout<<"DTOFPoint_factory: loading values from TOF data base"<<endl;
  } else {
    cout << "DTOFPoint_factory: Error loading values from TOF data base" <<endl;

    VELOCITY = 15.;  // set to some reasonable value
    HALFPADDLE = 126;   // set to some reasonable value
    BARWIDTH = 6.;
    return NOERROR;
  }

  VELOCITY    =    tofparms["TOF_C_EFFECTIVE"];
  HALFPADDLE  =    tofparms["TOF_HALFPADDLE"];
  BARWIDTH    =    tofparms["TOF_PADDLEWIDTH"];

  return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t DTOFPoint_factory::evnt(JEventLoop *loop, int eventnumber)
{
  // Get TOF hits
  vector<const DTOFHit*> hits;
  loop->Get(hits);

  // Local vectors for separating hits by plane
  vector<const DTOFHit *>uhits;
  vector<const DTOFHit *>vhits;
  for (unsigned int j = 0; j < hits.size(); j++){
    const DTOFHit *hit=hits[j];

    if (hit->orientation) uhits.push_back(hit);
    else vhits.push_back(hit);
  }

  // Match hits in the two planes by bar number and position as determined from
  // the time difference between two ends, where available.
  for (unsigned int i=0;i<uhits.size();i++){
    double ux    = uhits[i]->timediff;
    double utof  = uhits[i]->meantime;
    double dy    = uhits[i]->bar > 20 ? BARWIDTH : -BARWIDTH;
    double uy    = BARWIDTH*(uhits[i]->bar-20.5)+dy;
    double x_cut = BARWIDTH*0.55; // changed to "*0.55" from "/2." 4/28/2010  DL
    double t_cut = 1.;
    
    // Handle the single-ended bars
    if (uhits[i]->bar>40){
      x_cut=120.;
      t_cut=10.;
      switch(uhits[i]->bar){
      case 41:
	uy=-3.;
	ux=-72.;
	utof=uhits[i]->t_south-HALFPADDLE/VELOCITY;
	break;
      case 43:
	uy=3.;
	ux=-72.;
	utof=uhits[i]->t_south-HALFPADDLE/VELOCITY;
	break;
      case 42:
	uy=-3.;	
	ux=+72.;
	utof=uhits[i]->t_north-HALFPADDLE/VELOCITY;
	break;
      case 44:
	uy=+3.;
	ux=+72.;
	utof=uhits[i]->t_north-HALFPADDLE/VELOCITY;
	break;
      default:
	break;
      }

    }

    for (unsigned int j=0;j<vhits.size();j++){
      double vy    = vhits[j]->timediff;
      double vtof  = vhits[j]->meantime;
      double dx    = vhits[j]->bar > 20 ?  BARWIDTH:-BARWIDTH;
      double vx    = BARWIDTH*(vhits[j]->bar-20.5)+dx;
      double y_cut = BARWIDTH*0.55; // changed to "*0.55" from "/2." 4/28/2010  DL

      // Handle the single-ended bars
      if (vhits[j]->bar>40){
	y_cut=120.;
	t_cut=10.;
	switch(vhits[j]->bar){
	case 41:
	  vx=-3.;
	  vy=72.;
	  vtof=vhits[j]->t_south-HALFPADDLE/VELOCITY;
	  break;
	case 43:
	  vx=3.;
	  vy=72.;
	  vtof=vhits[j]->t_south-HALFPADDLE/VELOCITY;
	  break;
	case 42:
	  vx=-3.;
	  vy=-72.;
	  vtof=vhits[j]->t_north-HALFPADDLE/VELOCITY;
	  break;
	case 44:
	  vx=3.;
	  vy=-72;
	  vtof=vhits[j]->t_north-HALFPADDLE/VELOCITY;
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
	   point->pos.SetXYZ(ux,vy,618.81);
	   point->t=(utof+vtof)/2.;
	   point->dedx=(sqrt(uhits[i]->E_north*uhits[i]->E_south)
			+sqrt(vhits[j]->E_north*vhits[j]->E_south))/2.;
	 }
	 else if (vhits[j]->bar<41 && uhits[i]->bar>40){   
	   point->pos.SetXYZ(vx,vy,617.52);
	   point->t=vtof;
	   point->dedx=sqrt(vhits[j]->E_north*vhits[j]->E_south);
	 }
	 else{
	   point->pos.SetXYZ(ux,uy,620.10);
	   point->t=utof;
	   point->dedx=sqrt(uhits[i]->E_north*uhits[i]->E_south);
	 }
      	 _data.push_back(point);
      }
    }
  }
  
  return NOERROR;
}

