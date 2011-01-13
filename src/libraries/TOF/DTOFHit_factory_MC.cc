//    File: DTOFHit_factory_MC.cc
// Created: Thu Aug  9 11:56:15 EDT 2007
// Creator: B. Zihlmann
// comment: DTOFHit_factory_MC.cc replaces DTOFMCHit_factory.cc using the TAG "MC"
//          the DTOFHit.h data structure will be the same for real data and MC data
//

#include "DTOFHit_factory_MC.h"
#include "DTOFMCResponse.h"

#include <math.h>
using namespace std;

#define NaN std::numeric_limits<double>::quiet_NaN()

//------------------
// evnt
//------------------
jerror_t DTOFHit_factory_MC::evnt(JEventLoop *eventLoop, int eventnumber)
{

  vector<const DTOFMCResponse*> mcresponses;
  eventLoop->Get(mcresponses);

  for (unsigned int i = 0; i < mcresponses.size(); i++){

    const DTOFMCResponse *mcresponse = mcresponses[i];
    DTOFHit *hit = new DTOFHit;

    // keep basic information
    hit->id          = mcresponse->id;
    hit->orientation = mcresponse->orientation;
    hit->bar         = mcresponse->bar;
    
    // calculate meantime and time difference of MC data
    if (mcresponse->TDC_north>0 && mcresponse->TDC_south>0){
      hit->id          = mcresponse->id;
      hit->orientation = mcresponse->orientation;
      hit->bar = mcresponse->bar;

      float tn =  float(mcresponse->TDC_north)*TDC_RES_MC;
      float ts =  float(mcresponse->TDC_south)*TDC_RES_MC;
      float en =  float(mcresponse->ADC_north);
      float es =  float(mcresponse->ADC_south);
      // mean time
      float tm = (tn+ts)/2.;
      // time difference south-north so positive values are hits closer to north
      float td = (ts-tn);

      // position 
      float pos = C_EFFECTIVE*td/2.;
      
      hit->t_north = tn;
      hit->t_south = ts;

      hit->meantime  = tm;
      hit->timediff  = td;
      hit->pos       = pos;
      hit->dpos      = TOF_POS_RES; // see above only true if hit seen o both sides

      // mean energy deposition at the location of the hit position
      // devide by two to be comparable with single PMT hits
      en *= exp((HALFPADDLE-pos)/ATTEN_LENGTH) ;
      es *= exp((HALFPADDLE+pos)/ATTEN_LENGTH) ;
      hit->E_north = en;
      hit->E_south = es;

      float emean = (en+es)/2.; 
                                                    
      if (emean>2048) emean = 2048; // overflow
      emean = emean*TOF_ADC_TO_E;
      hit->dE = emean;

    } else {

      float tn =  float(mcresponse->TDC_north)*TDC_RES_MC;
      float ts =  float(mcresponse->TDC_south)*TDC_RES_MC;
      float en =  float(mcresponse->ADC_north);
      float es =  float(mcresponse->ADC_south);

      hit->t_north = tn;
      hit->E_north = en;
      hit->t_south = ts;
      hit->E_south = es;
      hit->meantime = NaN;
      hit->timediff = NaN;
      hit->pos      = NaN;
      hit->dpos     = NaN;
      hit->dE       = NaN;	
    }


    _data.push_back(hit);

  }

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTOFHit_factory_MC::brun(JEventLoop *loop, int eventnumber)
{
  map<string, double> tofparms;
  
  if ( !loop->GetCalib("TOF/tof_parms", tofparms)){
    cout<<"DTOFHit_factory_MC: loading values from TOF data base"<<endl;
  } else {
    cout << "DTOFHit_factory_MC: Error loading values from TOF data base" <<endl;
  }
  
  ATTEN_LENGTH   =    tofparms["TOF_ATTEN_LENGTH"]; 
  C_EFFECTIVE    =    tofparms["TOF_C_EFFECTIVE"];
  TDC_RES_MC     =    tofparms["TOF_TDC_RES_MC"];
  HALFPADDLE     =    tofparms["TOF_HALFPADDLE"];
  TOF_POS_RES    =    tofparms["TOF_POS_RES"];
  TOF_ADC_TO_E   =    tofparms["TOF_ADC_TO_E"];
  
  return NOERROR;
  
}
