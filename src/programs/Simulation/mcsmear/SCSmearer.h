// Smearing class for start counter (SC)

#ifndef _SCSMEARER_H_
#define _SCSMEARER_H_

#include "Smearer.h"


class sc_config_t 
{
  public:
	sc_config_t(JEventLoop *loop);

	double START_SIGMA;
	double START_PHOTONS_PERMEV;
	double START_PADDLE_THRESHOLD;

};


class SCSmearer : public Smearer
{
  public:
	SCSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		sc_config = new sc_config_t(loop);
	}
	~SCSmearer() {
		delete sc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	sc_config_t  *sc_config;
};



#endif // _SCSMEARER_H_