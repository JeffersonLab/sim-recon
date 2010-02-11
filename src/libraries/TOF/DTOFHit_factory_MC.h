//    File: DTOFHit_factory_MC.h
// Created: Thu Aug  9 11:56:15 EDT 2007
// Creator: B.Zihlmann 
// comment: replacement for DTOFMCHit_factory.h: now DTOFHIT.h will
//          be the same for real data and MC data
//  

#ifndef _DTOFHit_factory_MC
#define _DTOFHit_factory_MC

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTOFHit.h"

class DTOFHit_factory_MC:public JFactory<DTOFHit>{
 public:
  DTOFHit_factory_MC(){};
  ~DTOFHit_factory_MC(){};
  const char* Tag(void){return "MC";} // Monte Carlo TAG
  
 private:
  double TDC_RES_MC ;      // TDC resolution in [ns]
  double C_EFFECTIVE ;     // effective signal speed in paddle same as in hitFTOF.c
  double ATTEN_LENGTH ;    // effective attenuation legth in paddle same as in hitFTOF.c
  double TOF_POS_RES  ;    // TOF position resolution [cm] has to be determined experimentally
  // but depends on timing resolution and might be position dependent
  
  double TOF_ADC_TO_E ;    // convert ADC to energy deposition in paddle needs to be detemined
  double HALFPADDLE;       // half the detector paddle length
  
 protected:
  jerror_t brun(JEventLoop *eventLoop, int eventnumber);  ///< Called every event.
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);  ///< Called every event.
};

#endif // _DTOFHit_factory_MC

