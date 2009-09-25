// $Id$
//
//    File: DTOFMCResponse_factory.h
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DTOFMCResponse_factory_
#define _DTOFMCResponse_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTOFMCResponse.h"


class DTOFMCResponse_factory:public JFactory<DTOFMCResponse>{
 public:
  DTOFMCResponse_factory(){};
  ~DTOFMCResponse_factory(){};

 private:
  double ATTEN_LENGTH  ;   // attenuation length in paddle
  double C_EFFECTIVE   ;   // effective light speed in paddle
  double TWO_HIT_RESOL ;   // two hit timing resolution in paddle
  double THRESH_MEV    ;   // energy threshold in paddle
  double PHOTONS_PERMEV;   // scintillation photon generated per MeV energy deposition 
  double THETA_MAX     ;   // total internal reflection
  double PMT_SURFACE   ;   // PMT surface in cm^2
  double REFLECT_EFF   ;   // reflection efficiency of light 
  double PHE_EFF       ;   // efficiency to create photo elecrons
  double PMT_GAIN      ;   // PMT gain factor 4*10^7
  double ECHARGE       ;   // electric charge in C
  double ADC_RES       ;   // adc resolution pC/count
  double TOF_CENT_TRES ;   // time resolution of tof paddel for hit in the center log10(0.2ns)
  double TDC_RES       ;   // TDC resolution in ns
  double HALFPADDLE    ;   // half length of paddle
  
 private:
  
  jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DTOFMCResponse_factory_

