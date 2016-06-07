// main class for smearing data

#ifndef _SMEAR_H_
#define _SMEAR_H_

#include "HDDM/hddm_s.hpp"

#include <JANA/JEventLoop.h>
using namespace jana;

#include "mcsmear_config.h"
#include "smear_bcal.h"
using namespace bcal_smearing;


class Smear
{
    public:
		Smear(mcsmear_config_t *in_config, JEventLoop *loop) {   // constructor
			// create configuration classes
			config = in_config;
			bcal_config = new bcal_config_t(loop);
			fcal_config = new fcal_config_t(loop);
			cdc_config  = new cdc_config_t(loop);
			fdc_config  = new fdc_config_t(loop);
			ftof_config = new ftof_config_t(loop);
			sc_config   = new sc_config_t(loop);
			tagh_config = new tagh_config_t(loop);
			tagm_config = new tagm_config_t(loop);
			ps_config   = new ps_config_t(loop);
			psc_config  = new psc_config_t(loop);
			fmwpc_config = new fmwpc_config_t(loop);
			ccal_config = new ccal_config_t(loop);
		
			dBCALSmearer = new BCALSmearer(config, bcal_config);
		}
		~Smear() {  // destructor 
			delete dBCALSmearer;
		
			delete bcal_config;
			delete fcal_config;
			delete cdc_config;
			delete fdc_config;
			delete ftof_config;
			delete sc_config;
			delete tagh_config;
			delete tagm_config;
			delete ps_config;
			delete psc_config;
			delete fmwpc_config;
			delete ccal_config;
		}

		// main entrance - takes an event and smears it
		void SmearEvent(hddm_s::HDDM *record);

    private:
    	// utility functions
		void SetSeeds(const char *vals);
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
		bcal_config_t  *bcal_config;
		fcal_config_t  *fcal_config;
		cdc_config_t   *cdc_config;
		fdc_config_t   *fdc_config;
		ftof_config_t  *ftof_config;
		sc_config_t    *sc_config;
		tagh_config_t  *tagh_config;
		tagm_config_t  *tagm_config;
		
		ps_config_t    *ps_config;
		psc_config_t   *psc_config;
		fmwpc_config_t *fmwpc_config;
		ccal_config_t  *ccal_config;
		
};

#endif  // _SMEAR_H_