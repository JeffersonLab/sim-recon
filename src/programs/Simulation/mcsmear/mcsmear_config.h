// Classes to store configuration information for mcsmear

#ifndef _MCSMEAR_CONFIG_H_
#define _MCSMEAR_CONFIG_H_

#include "units.h"

#include <DANA/DApplication.h>
#include <JANA/JEventLoop.h>
#include "DRandom2.h"

using namespace jana;


// external function definitions from SampleGaussian.cc
double SampleGaussian(double sigma);
double SamplePoisson(double lambda);
double SampleRange(double x1, double x2);


// Overall configuration parameters
class mcsmear_config_t 
{
  public:
	mcsmear_config_t();
	~mcsmear_config_t();

	//-----------
	// SetSeeds
	//-----------
	void SetSeeds(const char *vals);

	// member variables
	bool ADD_NOISE;
	bool DROP_TRUTH_HITS;
	bool SMEAR_HITS;
    bool DUMP_RCDB_CONFIG;

	//bool SMEAR_BCAL;
	//bool FDC_ELOSS_OFF;
	bool IGNORE_SEEDS;
	double TRIGGER_LOOKBACK_TIME;
	bool APPLY_EFFICIENCY_CORRECTIONS;
	bool APPLY_HITS_TRUNCATION;

    bool FCAL_ADD_LIGHTGUIDE_HITS;
	
	// flags to pass command line info to subdetector classes
	double BCAL_NO_T_SMEAR;
	double BCAL_NO_DARK_PULSES;
	double BCAL_NO_SAMPLING_FLUCTUATIONS;
	double BCAL_NO_SAMPLING_FLOOR_TERM;
	double BCAL_NO_POISSON_STATISTICS;
	double BCAL_NO_FADC_SATURATION;
	double BCAL_NO_SIPM_SATURATION;
	
	
#ifdef HAVE_RCDB
    void LoadRCDBConnection();
	bool ParseRCDBConfigFile(int runNumber);
#endif  // HAVE_RCDB

    std::map<std::string, std::map<std::string, double> > readout;
};


#endif  // _MCSMEAR_CONFIG_H_
