// Smearing class for forward time-of-flight wall (TOF)

#ifndef _TOFSMEARER_H_
#define _TOFSMEARER_H_

#include "Smearer.h"

#include <TMath.h>

class tof_config_t 
{
  public:
	tof_config_t(JEventLoop *loop);
	
	inline double GetPaddleTimeResolution(int plane, int bar)  { 
		int paddle = plane*44 + bar - 1;   // hardcode for now
		return TOF_PADDLE_TIME_RESOLUTIONS.at(paddle); 
	}
	inline double GetHitTimeResolution(int plane, int bar)  { 
		// assume that the paddle resolution is given by: paddle resol = (hit resol)^2
		return GetPaddleTimeResolution(plane, bar)/TMath::Sqrt2(); 
	}

	double TOF_SIGMA;
	double TOF_PHOTONS_PERMEV;
	double TOF_BAR_THRESHOLD;

  	vector<double> TOF_PADDLE_TIME_RESOLUTIONS;
	
	vector< vector< pair<double,double> > > channel_efficiencies;

	double GetEfficiencyCorrectionFactor(hddm_s::FtofTruthHitList::iterator &siter) {
		if(siter->getEnd() == 0)
			return channel_efficiencies.at(siter->getPlane()).at(siter->getBar()-1).first;
		else 
			return channel_efficiencies.at(siter->getPlane()).at(siter->getBar()-1).second;
	}
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
