// Smearing class for tagger hodoscope (TAGH)

#ifndef _TAGHSMEARER_H_
#define _TAGHSMEARER_H_

#include "Smearer.h"


class tagh_config_t 
{
  public:
	tagh_config_t(JEventLoop *loop);

	double TAGH_TSIGMA;
	double TAGH_FADC_TSIGMA;
	double TAGH_NPE_PER_GEV;

};


class TAGHSmearer : public Smearer
{
  public:
	TAGHSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		tagh_config = new tagh_config_t(loop);
	}
	~TAGHSmearer() {
		delete tagh_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	tagh_config_t  *tagh_config;
};




#endif // _TAGHSMEARER_H_