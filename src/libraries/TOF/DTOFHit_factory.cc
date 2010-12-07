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

  return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t DTOFHit_factory::evnt(JEventLoop *loop, int eventnumber)
{

  vector<const DTOFHitRaw*> mcresponses;
  loop->Get(mcresponses);

  for (unsigned int i = 0; i < mcresponses.size(); i++){

    const DTOFHitRaw *mcresponse = mcresponses[i];
    DTOFHit *hit = new DTOFHit;

    // do any run dependent calibration here

    hit->id          = mcresponse->id;
    hit->orientation = mcresponse->plane;
    hit->t_north     = mcresponse->t_north;
    hit->E_north     = mcresponse->dE_north;
    hit->t_south     = mcresponse->t_south;
    hit->E_south     = mcresponse->dE_south;

    hit->meantime = (hit->t_north+hit->t_south)/2. - HALFPADDLE/C_EFFECTIVE;
    hit->timediff = (hit->t_south - hit->t_north)/2.*C_EFFECTIVE;

    _data.push_back(hit);

  }

  return NOERROR;
}

