// Smearing class for coarse pair spectrometer counters (PSC)

#ifndef _PSCSMEARER_H_
#define _PSCSMEARER_H_

#include "Smearer.h"


class psc_config_t 
{
  public:
	psc_config_t(JEventLoop *loop);

	double PSC_SIGMA;
	double PSC_PHOTONS_PERMEV;
	double PSC_THRESHOLD;

};


class PSCSmearer : public Smearer
{
  public:
	PSCSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		psc_config = new psc_config_t(loop);
	}
	~PSCSmearer() {
		delete psc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	psc_config_t  *psc_config;
};



#endif // _PSSMEARER_H_