// $Id$
//
//    File: DTOFPaddleHit_factory.cc
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//
// Modified: Wed Feb 12 13:19:10 EST 2014 by B. Zihlmann
//           reflect the changes in the TOF geometry with
//           19 LWB, 2 LNB, 2 SB, 2LNB, 19 LWB
//                          2 SB
//           LWB: long wide bars
//           LNB: long narrow bars
//           SB:  short bars
//           the bar numbering goes from 1 all through 46 with
//           bar 22 and 23 are the 4 short bars distinguished by north/south
//

#include <iostream>
using namespace std;

#include "DTOFPaddleHit_factory.h"
#include "DTOFHit.h"
#include "DTOFHitMC.h"
#include "DTOFPaddleHit.h"
#include <math.h>

//#define NaN std::numeric_limits<double>::quiet_NaN()
#define NaN __builtin_nan("")

//------------------
// brun
//------------------
jerror_t DTOFPaddleHit_factory::brun(JEventLoop *loop, int runnumber)
{

  map<string, double> tofparms;
 
  if ( !loop->GetCalib("TOF/tof_parms", tofparms)){
    cout<<"DTOFPaddleHit_factory: loading values from TOF data base"<<endl;
  } else {
    cout << "DTOFPaddleHit_factory: Error loading values from TOF data base" <<endl;

    C_EFFECTIVE = 15.;    // set to some reasonable value
    HALFPADDLE = 126;     // set to some reasonable value
    E_THRESHOLD = 0.0005; // energy threshold in GeV
    ATTEN_LENGTH = 400.;  // 400cm attenuation length
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
jerror_t DTOFPaddleHit_factory::evnt(JEventLoop *loop, int eventnumber)
{

  vector<const DTOFHit*> hits;
  loop->Get(hits,TOF_POINT_TAG.c_str());

  vector<const DTOFHit*> P1hitsL;
  vector<const DTOFHit*> P1hitsR;
  vector<const DTOFHit*> P2hitsL;
  vector<const DTOFHit*> P2hitsR;

  //int P1L[100];
  //int P1R[100];
  //int P2L[100];
  //int P2R[100];

  //int c1l = 0;
  //int c1r = 0;
  //int c2l = 0;
  //int c2r = 0;

  // sort the tof hits into left and right PMTs for both planes

  for (unsigned int i = 0; i < hits.size(); i++){
    const DTOFHit *hit = hits[i];
    if (hit->plane){
      if (hit->end){
	P2hitsR.push_back(hit);
	//P2R[c2r++] = i;
      } else {
	P2hitsL.push_back(hit);	
	//P2L[c2l++] = i;
      }
    } else {
      if (hit->end){
	P1hitsR.push_back(hit);
	//P1R[c1r++] = i;
      } else {
	P1hitsL.push_back(hit);
	//P1L[c1l++] = i;
      }
    }
  }

  for (unsigned int i=0; i<P1hitsL.size(); i++){
    int bar = P1hitsL[i]->bar;
    DTOFPaddleHit *hit = new DTOFPaddleHit;
    hit->bar = bar;
    hit->orientation   = P1hitsL[i]->plane;
    hit->E_north = P1hitsL[i]->dE;
    hit->t_north = P1hitsL[i]->t;
    hit->AddAssociatedObject(P1hitsL[i]);
    hit->E_south = 0.;
    hit->t_south = 0.;      
 
    if ((bar < 22) || (bar > 23)) {
      for (unsigned int j=0; j<P1hitsR.size(); j++){      
	if (bar==P1hitsR[j]->bar){
	  hit->E_south = P1hitsR[j]->dE;
	  hit->t_south = P1hitsR[j]->t;      
	  hit->AddAssociatedObject(P1hitsR[j]);
	}
      }
    }
    _data.push_back(hit);
  }


  for (unsigned int i=0; i<P1hitsR.size(); i++){   
    int bar = P1hitsR[i]->bar;
    int found = 0;

    if ((bar < 22) || (bar > 23)) {
      for (unsigned int j=0; j<P1hitsL.size(); j++){      
	if (bar==P1hitsL[j]->bar){
	  found = 1;
	}
      }
    }

    if (!found){
      DTOFPaddleHit *hit = new DTOFPaddleHit;
      hit->bar = bar;
      hit->orientation   = P1hitsR[i]->plane;
      hit->E_south = P1hitsR[i]->dE;
      hit->t_south = P1hitsR[i]->t;
      hit->E_north = 0.;
      hit->t_north = 0.;      
      hit->AddAssociatedObject(P1hitsR[i]);
      _data.push_back(hit);
    }
  }

  for (unsigned int i=0; i<P2hitsL.size(); i++){
    int bar = P2hitsL[i]->bar;
    DTOFPaddleHit *hit = new DTOFPaddleHit;
    hit->bar = bar;
    hit->orientation   = P2hitsL[i]->plane;
    hit->E_north = P2hitsL[i]->dE;
    hit->t_north = P2hitsL[i]->t;
    hit->AddAssociatedObject(P2hitsL[i]);
    hit->E_south = 0.;
    hit->t_south = 0.;      
 
    if ((bar < 22) || (bar > 23)) {
      for (unsigned int j=0; j<P2hitsR.size(); j++){      
	if (bar==P2hitsR[j]->bar){
	  hit->E_south = P2hitsR[j]->dE;
	  hit->t_south = P2hitsR[j]->t;      
	  hit->AddAssociatedObject(P2hitsR[j]);
	}
      }
    }
    _data.push_back(hit);
  }


  for (unsigned int i=0; i<P2hitsR.size(); i++){   
    int bar = P2hitsR[i]->bar;
    int found = 0;

    if ((bar < 22) || (bar > 23)) {
      for (unsigned int j=0; j<P2hitsL.size(); j++){      
	if (bar==P2hitsL[j]->bar){
	  found = 1;
	}
      }
    }

    if (!found){
      DTOFPaddleHit *hit = new DTOFPaddleHit;
      hit->bar = bar;
      hit->orientation   = P2hitsR[i]->plane;
      hit->E_south = P2hitsR[i]->dE;
      hit->t_south = P2hitsR[i]->t;
      hit->E_north = 0.;
      hit->t_north = 0.;      
      hit->AddAssociatedObject(P2hitsR[i]);
      _data.push_back(hit);
    }
  }


  for (int i=0;i<(int)_data.size(); i++) {
    
    DTOFPaddleHit *hit = _data[i];

    int check = -1;
    if (hit->E_north > E_THRESHOLD) {
      check++;
    }
    if (hit->E_south > E_THRESHOLD) {
      check++;
    }
    
    if (check > 0 ){
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
      hit->meantime = NaN;
      hit->timediff = NaN;
      hit->pos = NaN;
      hit->dpos = NaN;
      hit->dE = NaN;
   }

  }
  
  return NOERROR;
}

