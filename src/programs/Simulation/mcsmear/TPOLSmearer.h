// Smearing class for fine pair spectrometer counters (TPOL)

#ifndef _TPOLSMEARER_H_
#define _TPOLSMEARER_H_

#include "Smearer.h"


class tpol_config_t 
{
  public:
	tpol_config_t(JEventLoop *loop);

	double TPOL_SIGMA_NS;
	double TPOL_SIGMA1_MEV;
        double TPOL_SIGMA2_MEV;
	double TPOL_THRESHOLD_MEV;

};


class TPOLSmearer : public Smearer
{
  public:
	TPOLSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		tpol_config = new tpol_config_t(loop);
	}
	~TPOLSmearer() {
		delete tpol_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	tpol_config_t  *tpol_config;
};


#endif // _TPOLSMEARER_H_
