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
jerror_t DTOFPaddleHit_factory::brun(JEventLoop *loop, int32_t runnumber)
{

  map<string, double> tofparms;
 
  if ( !loop->GetCalib("TOF/tof_parms", tofparms)){
    //cout<<"DTOFPaddleHit_factory: loading values from TOF data base"<<endl;

    C_EFFECTIVE    =    tofparms["TOF_C_EFFECTIVE"];
    HALFPADDLE     =    tofparms["TOF_HALFPADDLE"];
    E_THRESHOLD    =    tofparms["TOF_E_THRESHOLD"];
    ATTEN_LENGTH   =    tofparms["TOF_ATTEN_LENGTH"];
  } else {
    cout << "DTOFPaddleHit_factory: Error loading values from TOF data base" <<endl;

    C_EFFECTIVE = 15.;    // set to some reasonable value
    HALFPADDLE = 126;     // set to some reasonable value
    E_THRESHOLD = 0.0005; // energy threshold in GeV
    ATTEN_LENGTH = 400.;  // 400cm attenuation length
  }

  ENERGY_ATTEN_FACTOR=exp(HALFPADDLE/ATTEN_LENGTH);
  TIME_COINCIDENCE_CUT=2.*HALFPADDLE/C_EFFECTIVE;

  if(loop->GetCalib("TOF/propagation_speed", propagation_speed))
    jout << "Error loading /TOF/propagation_speed !" << endl;

  if (loop->GetCalib("TOF/attenuation_lengths",AttenuationLengths))
    jout << "Error loading /TOF/attenuation_lengths !" <<endl;


  loop->Get(TOFGeom);

  return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t DTOFPaddleHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
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
    if (hit->has_fADC && hit->has_TDC){ // good hits have both ADC and TDC info
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
  }

  for (unsigned int i=0; i<P1hitsL.size(); i++){
    int bar = P1hitsL[i]->bar;
    if ((bar < TOFGeom[0]->FirstShortBar ) || (bar > TOFGeom[0]->LastShortBar)) {
      for (unsigned int j=0; j<P1hitsR.size(); j++){      
	if (bar==P1hitsR[j]->bar 
	    && fabs(P1hitsR[j]->t-P1hitsL[i]->t)<TIME_COINCIDENCE_CUT
	    && (P1hitsL[i]->dE>E_THRESHOLD || P1hitsR[j]->dE>E_THRESHOLD)){
	  DTOFPaddleHit *hit = new DTOFPaddleHit;
	  hit->bar = bar;
	  hit->orientation   = P1hitsL[i]->plane;
	  hit->E_north = P1hitsL[i]->dE;
	  hit->t_north = P1hitsL[i]->t;
	  hit->AddAssociatedObject(P1hitsL[i]);
	  hit->E_south = P1hitsR[j]->dE;
	  hit->t_south = P1hitsR[j]->t;      
	  hit->AddAssociatedObject(P1hitsR[j]);  

	  _data.push_back(hit);
	}
      }
    } 
  }
  
  for (unsigned int i=0; i<P1hitsL.size(); i++){ 
      int bar = P1hitsL[i]->bar;
      int found = 0;
      
      if ((bar < TOFGeom[0]->FirstShortBar) || (bar > TOFGeom[0]->LastShortBar)) {
      for (unsigned int j=0; j<P1hitsR.size(); j++){      
	if (bar==P1hitsR[j]->bar){
	  found = 1;
	}
      }
    }

    if (!found){
      if (P1hitsL[i]->dE>E_THRESHOLD){
	DTOFPaddleHit *hit = new DTOFPaddleHit;
	hit->bar = bar;
	hit->orientation   = P1hitsL[i]->plane;
	hit->E_north = P1hitsL[i]->dE;
	hit->t_north = P1hitsL[i]->t;
	hit->E_south = 0.;
	hit->t_south = 0.;  
	hit->AddAssociatedObject(P1hitsL[i]);

	_data.push_back(hit);
      }
    }
  }


	 for (unsigned int i=0; i<P1hitsR.size(); i++){   
	   int bar = P1hitsR[i]->bar;
	   int found = 0;
	   
	   if ((bar < TOFGeom[0]->FirstShortBar) || (bar > TOFGeom[0]->LastShortBar)) {
	     for (unsigned int j=0; j<P1hitsL.size(); j++){      
	       if (bar==P1hitsL[j]->bar){
		 found = 1;
	}
      }
    }

    if (!found){
      if (P1hitsR[i]->dE>E_THRESHOLD){
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
  }

  for (unsigned int i=0; i<P2hitsL.size(); i++){
    int bar = P2hitsL[i]->bar; 
    if ((bar <  TOFGeom[0]->FirstShortBar) || (bar > TOFGeom[0]->LastShortBar )){
      for (unsigned int j=0; j<P2hitsR.size(); j++){      
	if (bar==P2hitsR[j]->bar 
	    && fabs(P2hitsR[j]->t-P2hitsL[i]->t)<TIME_COINCIDENCE_CUT
	    && (P2hitsL[i]->dE>E_THRESHOLD || P2hitsR[j]->dE>E_THRESHOLD)){
	  DTOFPaddleHit *hit = new DTOFPaddleHit;
	  hit->bar = bar;
	  hit->orientation   = P2hitsL[i]->plane;
	  hit->E_north = P2hitsL[i]->dE;
	  hit->t_north = P2hitsL[i]->t;
	  hit->AddAssociatedObject(P2hitsL[i]);
	  hit->E_south = P2hitsR[j]->dE;
	  hit->t_south = P2hitsR[j]->t;      
	  hit->AddAssociatedObject(P2hitsR[j]);
	  
	  _data.push_back(hit);
	}
      }
    }
  }
  
	 for (unsigned int i=0; i<P2hitsL.size(); i++){   
	   int bar = P2hitsL[i]->bar;
	   int found = 0;

    if ((bar < TOFGeom[0]->FirstShortBar) || (bar > TOFGeom[0]->LastShortBar)) {
      for (unsigned int j=0; j<P2hitsR.size(); j++){      
	if (bar==P2hitsR[j]->bar){
	  found = 1;
	}
      }
    }

    if (!found){
      if (P2hitsL[i]->dE>E_THRESHOLD){
	DTOFPaddleHit *hit = new DTOFPaddleHit;
	hit->bar = bar;
	hit->orientation   = P2hitsL[i]->plane;
	hit->E_north = P2hitsL[i]->dE;
	hit->t_north = P2hitsL[i]->t;
	hit->E_south = 0.;
	hit->t_south = 0.;      
	hit->AddAssociatedObject(P2hitsL[i]);

	_data.push_back(hit);
      }
    }
  }


  for (unsigned int i=0; i<P2hitsR.size(); i++){   
    int bar = P2hitsR[i]->bar;
    int found = 0;

    if ((bar < TOFGeom[0]->FirstShortBar) || (bar > TOFGeom[0]->LastShortBar)) {
      for (unsigned int j=0; j<P2hitsL.size(); j++){      
	if (bar==P2hitsL[j]->bar){
	  found = 1;
	}
      }
    }

    if (!found){
       if (P2hitsR[i]->dE>E_THRESHOLD){
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
      int id=44*hit->orientation+hit->bar-1;
      double v=propagation_speed[id];
      hit->meantime = (hit->t_north+hit->t_south)/2. - HALFPADDLE/v;
      hit->timediff = (hit->t_south - hit->t_north)/2.;
      float pos = hit->timediff * v;  
      hit->pos = pos;
      hit->dpos      = 2.;  // manually/artificially set to 2cm. 
      
      // mean energy deposition at the location of the hit position
      // use geometrical mean
      //hit->dE = ENERGY_ATTEN_FACTOR*sqrt(hit->E_north*hit->E_south);

      float xl = 126. - pos; // distance to left PMT 
      float xr = 126. + pos; // distance to right PMT
      int idl = hit->orientation*88 + hit->bar-1;
      int idr = idl+44;
      float d1 = AttenuationLengths[idl][0];
      float d2 = AttenuationLengths[idl][1];
      // reference distance is 144cm from PMT
      float att_left = ( TMath::Exp(-144./d1) +  TMath::Exp(-144./d2)) / 
	( TMath::Exp(-xl/d1) +  TMath::Exp(-xl/d2));
      d1 = AttenuationLengths[idr][0];
      d2 = AttenuationLengths[idr][1];
      float att_right = ( TMath::Exp(-144./d1) +  TMath::Exp(-144./d2)) / 
	( TMath::Exp(-xr/d1) +  TMath::Exp(-xr/d2));
      hit->dE = (hit->E_north*att_left + hit->E_south*att_right)/2.;
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

