// Smearing class for central drift chamber (CDC)

#ifndef _CDCSMEARER_H_
#define _CDCSMEARER_H_

#include "Smearer.h"


class cdc_config_t 
{
  public:
	cdc_config_t(JEventLoop *loop);

	double CDC_TDRIFT_SIGMA;
	double CDC_TIME_WINDOW;
	double CDC_PEDESTAL_SIGMA;
	double CDC_THRESHOLD_FACTOR; // number of pedestal sigmas for determining sparsification threshold

	vector< vector<double> > wire_efficiencies;

	void CalcNstraws(JEventLoop *loop, int32_t runnumber, vector<unsigned int> &Nstraws);
	double GetEfficiencyCorrectionFactor(int ring, int straw) {
		return wire_efficiencies.at(ring-1).at(straw-1);
	}
};


class CDCSmearer : public Smearer
{
  public:
	CDCSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		cdc_config = new cdc_config_t(loop);
	}
	~CDCSmearer() {
		delete cdc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	cdc_config_t  *cdc_config;
};



#endif // _CDCSMEARER_H_
