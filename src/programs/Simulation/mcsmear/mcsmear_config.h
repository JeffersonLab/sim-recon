// Classes to store configuration information for mcsmear

#ifndef _MCSMEAR_CONFIG_H_
#define _MCSMEAR_CONFIG_H_

#include "units.h"

#include <DANA/DApplication.h>
#include <JANA/JEventLoop.h>
#include "DRandom2.h"

#include <RCDB/Connection.h>

using namespace jana;

//  std::numeric_limits::epsilon();

// external function definitions from SampleGaussian.cc
double SampleGaussian(double sigma);
double SamplePoisson(double lambda);
double SampleRange(double x1, double x2);


// Overall configuration parameters
class mcsmear_config_t 
{
  public:
	mcsmear_config_t() {
		// default values
		DROP_TRUTH_HITS = false;
		ADD_NOISE      = false;
		SMEAR_HITS     = true;
		//SMEAR_BCAL     = true;
		IGNORE_SEEDS   = false;
		
		TRIGGER_LOOKBACK_TIME = -100; // ns
		
#ifdef RCDB_HOME
		// RCDB configuration
		// first determine which database to connect to
		if( getenv("RCDB_CONNECTION")!= NULL )
			RCDB_CONNECTION = getenv("RCDB_CONNECTION");
		else
			RCDB_CONNECTION = "mysql://rcdb@hallddb/rcdb";   // default to outward-facing MySQL DB
			
		gPARMS->SetDefaultParameters("RCDB_CONNECTION", RCDB_CONNECTION, "URL used to access RCDB.");
			
		// load connection to RCDB
		rcdb_connection = new rcdb::Connection(RCDB_CONNECTION);
#endif //RCDB_HOME
	}
	~mcsmear_config_t() {}

	//-----------
	// SetSeeds
	//-----------
	void SetSeeds(const char *vals) {
   		/// This is called from the command line parser to
   		/// set the initial seeds based on user input from
   		/// the command line.
   		//
   		stringstream ss(vals);
   		Int_t seed1, seed2, seed3;
   		ss >> seed1 >> seed2 >> seed3;
   		UInt_t *useed1 = reinterpret_cast<UInt_t*>(&seed1);
   		UInt_t *useed2 = reinterpret_cast<UInt_t*>(&seed2);
   		UInt_t *useed3 = reinterpret_cast<UInt_t*>(&seed3);
   		gDRandom.SetSeeds(*useed1, *useed2, *useed3);

   		cout << "Seeds set from command line. Any random number" << endl;
   		cout << "seeds found in the input file will be ignored!" << endl;
   		IGNORE_SEEDS = true;
	}

	// member variables
	bool ADD_NOISE;
	bool DROP_TRUTH_HITS;
	bool SMEAR_HITS;
	//bool SMEAR_BCAL;
	//bool FDC_ELOSS_OFF;
	bool IGNORE_SEEDS;
	double TRIGGER_LOOKBACK_TIME;
	
#ifdef RCDB_HOME
	// RCDB information
	string RCDB_CONNECTION;
	rcdb::Connection *rcdb_connection;
#endif  // RCDB_HOME

};


#endif  // _MCSMEAR_CONFIG_H_
