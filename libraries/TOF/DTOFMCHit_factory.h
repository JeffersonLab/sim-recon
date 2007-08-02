// $Id: DTOFMCHit_factory.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DTOFMCHit_factory.h
// Created: Mon Jul  9 16:21:12 EDT 2007
// Creator: B.Zihlmann 
//

#ifndef _DTOFMCHit_factory_
#define _DTOFMCHit_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DTOFMCHit.h"

class DTOFMCHit_factory:public JFactory<DTOFMCHit>{
 public:
  DTOFMCHit_factory(){};
  ~DTOFMCHit_factory(){};
  const string toString(void);
  
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

#endif // _DTOFMCHit_factory_

