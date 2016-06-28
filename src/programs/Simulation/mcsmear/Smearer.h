// abstract base class for smearing hits in a subdetector

#ifndef _SMEARER_H_
#define _SMEARER_H_

#include "mcsmear_config.h"
#include "HDDM/hddm_s.hpp"

#include <JANA/JEventLoop.h>
using namespace jana;

class Smearer
{
  public:
	Smearer(JEventLoop *loop, mcsmear_config_t *in_config) {
		config = in_config;
	};
	
	virtual void SmearEvent(hddm_s::HDDM *record) = 0;
	
  protected:
  	mcsmear_config_t *config;  // save a link to this information, but we do not own it

};

#endif // _SMEARER_H_