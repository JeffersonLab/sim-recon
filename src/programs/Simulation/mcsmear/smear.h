// main class for smearing data

#ifndef _SMEAR_H_
#define _SMEAR_H_

#include "mcsmear_config.h"
#include "smear_bcal.h"

#include "HDDM/hddm_s.hpp"

class Smear
{
    public:
		Smear(mcsmear_config_t *in_config) {   // constructor
			// create configuration classes
			config = in_config;
		
			dBCALSmearer = new BCALSmearer(config, bcal_config);
		}
		~Smear() {}  // destructor 

		// main entrance - takes an event and smears it
		void Smear(hddm_s::HDDM *record);

    private:
    	// utility functions
		void SetSeeds(const char *vals)
		void GetAndSetSeeds(hddm_s::HDDM *record);

		// functions to smear individual detectors
		void SmearCDC(hddm_s::HDDM *record);
		void SmearFDC(hddm_s::HDDM *record);
		void SmearFCAL(hddm_s::HDDM *record);
		void SmearCCAL(hddm_s::HDDM *record);
		//void SmearBCAL(hddm_s::HDDM *record);
		void SmearTOF(hddm_s::HDDM *record);
		void SmearSTC(hddm_s::HDDM *record);
		void SmearCherenkov(hddm_s::HDDM *record);
		void SmearTAGM(hddm_s::HDDM *record);
		void SmearTAGH(hddm_s::HDDM *record);
		void SmearPS(hddm_s::HDDM *record);
		void SmearPSC(hddm_s::HDDM *record);
		void SmearFMWPC(hddm_s::HDDM *record);

		// More complicated smearing procedures are implemented in their own classes
		BCALSmearer *dBCALSmearer;

		// store configuration information
		mcsmear_config_t *config;
};

#endif  // _SMEAR_H_