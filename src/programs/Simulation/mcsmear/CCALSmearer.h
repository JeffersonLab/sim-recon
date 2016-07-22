// Smearing class for compton calorimeter (CCAL)

#ifndef _CCALSMEARER_H_
#define _CCALSMEARER_H_

#include "Smearer.h"

#include <CCAL/DCCALGeometry.h>


class ccal_config_t 
{
  public:
	ccal_config_t(JEventLoop *loop);

	// Time smearing factor
	double CCAL_SIGMA;

	// Photon-statistics factor for smearing hit energy for CompCal
	double CCAL_PHOT_STAT_COEF;

	// Single block energy threshold (applied after smearing)
	double CCAL_BLOCK_THRESHOLD;

};


class CCALSmearer : public Smearer
{
  public:
	CCALSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		ccal_config = new ccal_config_t(loop);
		ccalGeom = new DCCALGeometry();
	}
	~CCALSmearer() {
		delete ccal_config;
		delete ccalGeom;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	ccal_config_t  *ccal_config;
	DCCALGeometry  *ccalGeom;
};


#endif // _CCALSMEARER_H_