// Smearing class for tagger microscope (TAGM)

#ifndef _TAGMSMEARER_H_
#define _TAGMSMEARER_H_

#include "Smearer.h"


class tagm_config_t 
{
  public:
	tagm_config_t(JEventLoop *loop);

	double TAGM_TSIGMA;
	double TAGM_FADC_TSIGMA;
	double TAGM_NPIX_PER_GEV;

};


class TAGMSmearer : public Smearer
{
  public:
	TAGMSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		tagm_config = new tagm_config_t(loop);
	}
	~TAGMSmearer() {
		delete tagm_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	tagm_config_t  *tagm_config;
};



#endif // _TAGMSMEARER_H_