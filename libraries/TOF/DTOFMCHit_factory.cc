// $Id: DTOFMCHit_factory.cc 2710 2007-06-22 15:09:48Z zihlmann $
//
//    File: DTOFMCHit_factory.cc
// Created: Mon Jul  9 16:34:24 EDT 2007
// Creator: B. Zihlmann
//
#include "DTOFMCHit_factory.h"
#include "DTOFMCResponse.h"

#include <math.h>

#define TDC_MC_RES     0.06   // TDC resolution in [ns]
#define C_EFFECTIVE    15.0   // effective signal speed in paddle same as in hitFTOF.c
#define ATTEN_LENGTH   150    // effective attenuation legth in paddle same as in hitFTOF.c
#define TOF_POS_RES    3.0    // TOF position resolution [cm] has to be determined experimentally
                              // but depends on timing resolution and might be position dependent

#define TOF_ADC_TO_E      1.   // convert ADC to energy deposition in paddle needs to be detemined

// the following has to be done better by modifing DTOFGeometry class
#define HALFPADDLE     126.   // half paddle legth

void debugfuncMCHit(){
  //do nothing
}

//------------------
// evnt
//------------------
jerror_t DTOFMCHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{

  vector<const DTOFMCResponse*> mcresponses;
  eventLoop->Get(mcresponses);

  for (unsigned int i = 0; i < mcresponses.size(); i++){

    const DTOFMCResponse *mcresponse = mcresponses[i];
    DTOFMCHit *hit = new DTOFMCHit;

    // calculate meantime and time difference of MC data
    if (mcresponse->TDC_north>0 && mcresponse->TDC_south>0){
      hit->id          = mcresponse->id;
      hit->orientation = mcresponse->orientation;

      float tn =  float(mcresponse->TDC_north)*TDC_MC_RES;
      float ts =  float(mcresponse->TDC_south)*TDC_MC_RES;
      float en =  float(mcresponse->ADC_north);
      float es =  float(mcresponse->ADC_south);
      // energy weighted mean time
      float tm = (tn*en+ts*es)/(en+es)/2.;
      // time difference south-north so positive values are hits closer to north
      float td = (ts-tn);

      // position 
      float pos = C_EFFECTIVE*td/2.;
      
      hit->meantime  = tm;
      hit->timediff = td;
      hit->pos      = pos;
      hit->dpos     = TOF_POS_RES; // see above only true if hit seen o both sides

      // mean energy deposition weight by arrival time at PMT (favor closer PMT)
      // devide by two to be comparable with single PMT hits
      en *= exp((HALFPADDLE-pos)/ATTEN_LENGTH) ;
      es *= exp((HALFPADDLE+pos)/ATTEN_LENGTH) ;

      float emean = (en/tn+es/ts)/(1./tn+1./ts)/2.; 
                                                    
      emean = emean*TOF_ADC_TO_E;
      hit->dE = emean;

    } else {
      hit->id          = mcresponse->id;
      hit->orientation = mcresponse->orientation;
      hit->meantime = -999.;
      hit->timediff = -999.;
      hit->pos      = -999.;
      hit->dpos     = -999.;
      hit->dE       = -999.;	
    }

    debugfuncMCHit();

    _data.push_back(hit);

  }

  return NOERROR;
}


//------------------
// toString
//------------------
const string DTOFMCHit_factory::toString(void)
{
  // Ensure our Get method has been called so _data is up to date
  Get();
  if(_data.size()<=0)return string(); // don't print anything if we have no data!

  printheader( "id: orientation: pos[cm]:  epos[cm]:  dE [MeV]: meantime [ns]: timediff [ns]:" );

	
  for(unsigned int i=0; i<_data.size(); i++){
    DTOFMCHit *tofhit = _data[i];

    printnewrow();
    printcol("%d",	tofhit->id );
    printcol("%d",	tofhit->orientation );
    printcol("%2.3f",	tofhit->pos );
    printcol("%2.3f",	tofhit->dpos );
    printcol("%1.3f",	tofhit->dE );
    printcol("%1.3f",	tofhit->meantime );    
    printcol("%1.3f",	tofhit->timediff );

    printrow();
  }
  
	
  return _table;
}
