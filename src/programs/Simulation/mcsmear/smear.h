// main class for smearing data

#ifndef _SMEAR_H_
#define _SMEAR_H_

#include <map>
using namespace std;

#include "HDDM/hddm_s.hpp"
#include "GlueX.h"

#include <JANA/JEventLoop.h>
using namespace jana;

#include "mcsmear_config.h"

#include <Smearer.h>
#include <CDCSmearer.h>
#include <FDCSmearer.h>
#include <BCALSmearer.h>
#include <FCALSmearer.h>
#include <CCALSmearer.h>
#include <TOFSmearer.h>
#include <SCSmearer.h>
#include <FDIRCSmearer.h>
#include <TAGHSmearer.h>
#include <TAGMSmearer.h>
#include <PSSmearer.h>
#include <PSCSmearer.h>
#include <TPOLSmearer.h>
#include <FMWPCSmearer.h>



class Smear
{
    public:
		Smear(mcsmear_config_t *in_config, JEventLoop *loop, string detectors_to_load="all");
		~Smear();

		// main entrance - takes an event and smears it
		void SmearEvent(hddm_s::HDDM *record);

    private:
    	// utility functions
		void SetSeeds(const char *vals);
		void GetAndSetSeeds(hddm_s::HDDM *record);

		// Detector digitization/smearing is implemented in a different class for each subdetector
		map<DetectorSystem_t, Smearer *>  smearers;
		
		mcsmear_config_t *config;
};

#endif  // _SMEAR_H_
