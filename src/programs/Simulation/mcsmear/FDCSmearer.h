// Smearing class for forward drift chamber (FDC)

#ifndef _FDCSMEARER_H_
#define _FDCSMEARER_H_

#include "Smearer.h"


class fdc_config_t 
{
  public:
	fdc_config_t(JEventLoop *loop);

	double FDC_TDRIFT_SIGMA;
	double FDC_CATHODE_SIGMA;
	double FDC_PED_NOISE;         // in pC calculated in above
	double FDC_THRESHOLD_FACTOR;  // number of pedestal sigmas for determining sparcification threshold
	//double FDC_HIT_DROP_FRACTION = 0.0; // 1000.0E-9
	double FDC_TIME_WINDOW;
	double FDC_THRESH_KEV;        // fdc anode discriminator threshold

	vector< vector<double > > channel_efficiencies;
	
	double GetEfficiencyCorrectionFactor(hddm_s::FdcCathodeStripList::iterator &siter) {
		// cathode strips
		int gPlane = 9*(siter->getModule()-1) + 3*(siter->getLayer()-1) 
             		+ (siter->getPlane()-1);
        int element = siter->getStrip();

		return channel_efficiencies.at(gPlane-1).at(element-1);
	}

	double GetEfficiencyCorrectionFactor(hddm_s::FdcAnodeWireList::iterator &witer) {
		// anode wires
		int gPlane = 9*(witer->getModule()-1) + 3*(witer->getLayer()-1) + 1; 
        int element = witer->getWire();

		return channel_efficiencies.at(gPlane-1).at(element-1);
	}

};



class FDCSmearer : public Smearer
{
  public:
	FDCSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		fdc_config = new fdc_config_t(loop);
	}
	~FDCSmearer() {
		delete fdc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	fdc_config_t  *fdc_config;
};



#endif // _FDCSMEARER_H_