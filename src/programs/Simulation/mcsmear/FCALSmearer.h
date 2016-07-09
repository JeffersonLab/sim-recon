// Smearing class for forward calorimeter (FCAL)

#ifndef _FCALSMEARER_H_
#define _FCALSMEARER_H_

#include "Smearer.h"

#include <FCAL/DFCALGeometry.h>


class fcal_config_t 
{
  public:
	fcal_config_t(JEventLoop *loop, DFCALGeometry *fcalGeom);

	double FCAL_PHOT_STAT_COEF;
	double FCAL_BLOCK_THRESHOLD;
	double FCAL_TSIGMA;
	
	vector<double> FCAL_GAINS;
	double FCAL_MC_ESCALE;
	
	vector< vector<double > > block_efficiencies;
	
	double GetEfficiencyCorrectionFactor(double row, double column) {
		return block_efficiencies.at(row).at(column);
	}
};



class FCALSmearer : public Smearer
{
  public:
	FCALSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		fcalGeom = new DFCALGeometry();
		fcal_config = new fcal_config_t(loop, fcalGeom);
	}
	~FCALSmearer() {
		delete fcal_config;
		delete fcalGeom;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	fcal_config_t  *fcal_config;
  	DFCALGeometry *fcalGeom;
};


#endif // _FCALSMEARER_H_