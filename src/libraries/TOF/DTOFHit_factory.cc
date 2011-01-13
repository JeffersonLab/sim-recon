// $Id$
//
//    File: DTOFHit_factory.cc
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//
#include <iostream>
using namespace std;

#include "DTOFHit_factory.h"
#include "DTOFHitRaw.h"
#include "DTOFHit.h"
#include <math.h>


//------------------
// brun
//------------------
jerror_t DTOFHit_factory::brun(JEventLoop *loop, int runnumber)
{

  map<string, double> tofparms;
 
  if ( !loop->GetCalib("TOF/tof_parms", tofparms)){
    cout<<"DTOFHit_factory: loading values from TOF data base"<<endl;
  } else {
    cout << "DTOFHit_factory: Error loading values from TOF data base" <<endl;

    C_EFFECTIVE = 15.;  // set to some reasonable value
    HALFPADDLE = 126;   // set to some reasonable value
    return NOERROR;
  }

  C_EFFECTIVE    =    tofparms["TOF_C_EFFECTIVE"];
  HALFPADDLE     =    tofparms["TOF_HALFPADDLE"];
  E_THRESHOLD    =    tofparms["TOF_E_THRESHOLD"];
  ATTEN_LENGTH   =    tofparms["TOF_ATTEN_LENGTH"];
  return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t DTOFHit_factory::evnt(JEventLoop *loop, int eventnumber)
{

  vector<const DTOFHitRaw*> mcresponses;
  loop->Get(mcresponses,TOF_POINT_TAG.c_str());

  for (unsigned int i = 0; i < mcresponses.size(); i++){

    const DTOFHitRaw *mcresponse = mcresponses[i];
    DTOFHit *hit = new DTOFHit;

    // do any run dependent calibration here

    hit->id          = mcresponse->id;
    hit->orientation = mcresponse->plane;
    hit->bar         = mcresponse->bar;

    int check = -1;
    if (mcresponse->dE_north > E_THRESHOLD) {
      hit->t_north     = mcresponse->t_north;
      hit->E_north     = mcresponse->dE_north;
      check++;
    } else {
      hit->t_north  = -999.;
      hit->E_north  = -999.;
    }
    if (mcresponse->dE_south > E_THRESHOLD) {
      hit->t_south     = mcresponse->t_south;
      hit->E_south     = mcresponse->dE_south;
      check++;
    } else {
      hit->t_south     = -999.;
      hit->E_south     = -999.;
    }
    
    if (check > 0 ) {
      hit->meantime = (hit->t_north+hit->t_south)/2. - HALFPADDLE/C_EFFECTIVE;
      hit->timediff = (hit->t_south - hit->t_north)/2.;
      float pos = hit->timediff * C_EFFECTIVE;  
      hit->pos = pos;
      hit->dpos      = 2.;  // manually/artificially set to 2cm. 
  
    // mean energy deposition at the location of the hit position
      // devide by two to be comparable with single PMT hits
      float en = hit->E_north  * exp((HALFPADDLE-pos)/ATTEN_LENGTH) ;
      float es = hit->E_south  * exp((HALFPADDLE+pos)/ATTEN_LENGTH) ;
      float emean = (en+es)/2.; 
      hit->dE = emean;
 
    } else {
      hit->meantime = -999.;
      hit->timediff = -999.;
      hit->pos = -999.;
      hit->dpos = -999.;
      hit->dE = -999.;
   }

    _data.push_back(hit);

  }

  return NOERROR;
}

