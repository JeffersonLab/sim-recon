// $Id$
//
//    File: DMCTrigger_factory.cc
// Created: Tue Jun  7 10:15:05 EDT 2011  (originally was DTrigger_factory.cc)
// Creator: davidl (on Darwin eleanor.jlab.org 10.7.0 i386)
//


#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include <JANA/JApplication.h>
#include "DMCTrigger_factory.h"
using namespace jana;

#include <BCAL/DBCALHit.h>
#include <BCAL/DBCALDigiHit.h>
#include <FCAL/DFCALHit.h>
#include <START_COUNTER/DSCHit.h>
#include <TOF/DTOFHit.h>

// constant stolen from mcsmear - this should be put in the CCDB
static double BCAL_MEV_PER_ADC_COUNT    = 0.029;

//------------------
// init
//------------------
jerror_t DMCTrigger_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DMCTrigger_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
	// Get attenuation parameters
	double L_over_2 = DBCALGeometry::GetBCAL_length()/2.0;
	double Xo = DBCALGeometry::ATTEN_LENGTH;
	unattenuate_to_center = exp(+L_over_2/Xo);

    USE_OLD_BCAL_HITS = 0;
    gPARMS->SetDefaultParameter("TRIGGER:USE_OLD_BCAL_HITS", USE_OLD_BCAL_HITS,
                                "Use old DBCALHit calculation for BCAL energy sum.");

    // Per-channel trigger thresholds
    // Defaults in data according to A. Somov for Spring 2015 run (as of 11/13/15):
    // FCAL: 130 MeV (65 adc counts)
    // BCAL: 13.2 MeV (20 adc counts)
    // Added 11/17/15 (sdobbs):
    // Caveat: BCAL simulations use an adc count -> MeV conversion factor of 0.029, 
    // so for the BCAL use a threshold of  
    BCAL_CHANNEL_THRESHOLD = 455.;     // in adc counts
    FCAL_CHANNEL_THRESHOLD = 0.130;    // in MeV
    gPARMS->SetDefaultParameter("TRIGGER:BCAL_CHANNEL_THRESHOLD", BCAL_CHANNEL_THRESHOLD,
                                "Threshold for BCAL energy sum (units of adc counts)");
    gPARMS->SetDefaultParameter("TRIGGER:FCAL_CHANNEL_THRESHOLD", FCAL_CHANNEL_THRESHOLD,
                                "Threshold for FCAL energy sum (units of MeV)");


	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCTrigger_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// See comments in DMCTrigger_factory.h

    vector<const DBCALHit*> bcalhits;
	vector<const DFCALHit*> fcalhits;
	vector<const DSCHit*> schits;
	vector<const DTOFHit*> tofhits;
    loop->Get(bcalhits);
	loop->Get(fcalhits);
	loop->Get(schits);
	loop->Get(tofhits);

    double Ebcal = 0.0;
    double Ebcal_all = 0.0;
    double Efcal = 0.0;
    double Efcal_all = 0.0;
    double BCAL_Eupstream = 0.0;
    double BCAL_Edownstream = 0.0;
    // BCAL simulations now write out DBCALDigiHits, so we can directly sum the unattenuated energies
    // Keep the old algorithm around for compatabilty
    if(USE_OLD_BCAL_HITS) {
        /// In GlueX-doc-1043, it appaears the energy deposited in the BCAL
        /// is used. In reality, only the attenuated energy will be available
        /// to the L1 trigger electronics. If an average of the sum of both
        /// sides is used, it will be off by at most, 20% at the ends (excluding
        /// any dark pulsing effects). This conclusion comes from assuming
        /// we scale the values by e^(L/2/lambda) such that the value of the
        /// threshold is in GeV for energy deposited at the center of the module.
        /// The average of the two ends then is:
        ///
        ///   Eavg = (E*e^+x + E*e^-x)/2
        ///
        /// or
        ///
        ///   Eavg = E * cosh(x)
        ///
        /// where x = L/2/lambda 
        ///       L = 390 cm
        ///       lambda = 300 cm
        ///
        /// The dark pulses should contribute at most about 10 MeV.
        ///
        /// For the BCAL "energy" we therefore take the average of the
        /// sums of the DBCALHits for each side knowing it is
        /// overestimated by as much as 20%.
        ///
        /// The FCAL energy is just a straight sum of the FCAL hits
        
        for(unsigned int i=0; i< bcalhits.size(); i++){
            const DBCALHit* bcalhit = bcalhits[i];
            
            if(bcalhit->end == DBCALGeometry::kUpstream){
                Ebcal_all += bcalhit->E;
                if(1000.*bcalhit->E > BCAL_CHANNEL_THRESHOLD*BCAL_MEV_PER_ADC_COUNT) 
                    BCAL_Eupstream += bcalhit->E;
            }else{
                Ebcal_all += bcalhit->E;
                if(1000.*bcalhit->E > BCAL_CHANNEL_THRESHOLD*BCAL_MEV_PER_ADC_COUNT) 
                    BCAL_Edownstream += bcalhit->E;
            }
        }
	
        // Calculate "energy" sum for BCAL in GeV-ish units
        Ebcal = unattenuate_to_center * (BCAL_Eupstream + BCAL_Edownstream)/2.0;
	
        // Sum up FCAL energy
        for(unsigned int i=0; i< fcalhits.size(); i++){
            const DFCALHit* fcalhit = fcalhits[i];
	   
            Efcal_all += fcalhit->E;
            if(fcalhit->E > FCAL_CHANNEL_THRESHOLD)
                Efcal += fcalhit->E;
        }
    } else {
        // Sum up BCAL energy
        for(unsigned int i=0; i< bcalhits.size(); i++){
            const DBCALHit* bcalhit = bcalhits[i];
            const DBCALDigiHit* bcaldigihit = NULL;
            bcalhit->GetSingle(bcaldigihit);
            if(bcaldigihit == NULL)
                continue;
            
            Ebcal_all += bcalhit->E;
            if(bcaldigihit->pulse_integral > BCAL_CHANNEL_THRESHOLD)  
                Ebcal += bcalhit->E;
        }
        
        // Sum up FCAL energy
        for(unsigned int i=0; i< fcalhits.size(); i++){
            const DFCALHit* fcalhit = fcalhits[i];
	
            Efcal_all += fcalhit->E;
            if(fcalhit->E > FCAL_CHANNEL_THRESHOLD)
                Efcal += fcalhit->E;
        }
    }
        
	// Number of start counter hits
	unsigned int Nschits = schits.size();

	// Number of TOF hits
	unsigned int Ntofhits = tofhits.size();

	DMCTrigger *trig = new DMCTrigger;

	// BCAL and FCAL
	bool sum_cut = (Ebcal + 4.0*Efcal)>=2.0;
	trig->L1a_fired = sum_cut && Ebcal>0.200 && Efcal>0.030;
	trig->L1b_fired = sum_cut && Ebcal>0.030 && Efcal>0.030 && Nschits>0;
	trig->L1c_fired = Efcal>=0.250;

	trig->Ebcal = Ebcal;
	trig->Efcal = Efcal;
	trig->Ebcal_all = Ebcal_all;
	trig->Efcal_all = Efcal_all;
	trig->Nschits = Nschits;
	trig->Ntofhits = Ntofhits;
	
	_data.push_back(trig);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DMCTrigger_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DMCTrigger_factory::fini(void)
{
	return NOERROR;
}

