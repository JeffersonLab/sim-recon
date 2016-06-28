// Smearing class for forward time-of-flight wall (TOF)

#ifndef _TOFSMEARER_H_
#define _TOFSMEARER_H_

#include "Smearer.h"


class tof_config_t 
{
  public:
	tof_config_t(JEventLoop *loop);
	
	double TOF_SIGMA;
	double TOF_PHOTONS_PERMEV;
	double TOF_BAR_THRESHOLD;

};


class TOFSmearer : public Smearer
{
  public:
	TOFSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		tof_config = new tof_config_t(loop);
	}
	~TOFSmearer() {
		delete tof_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	tof_config_t  *tof_config;
};


#endif // _TOFSMEARER_H_