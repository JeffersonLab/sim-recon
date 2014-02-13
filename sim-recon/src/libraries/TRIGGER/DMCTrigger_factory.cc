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
#include <FCAL/DFCALHit.h>
#include <START_COUNTER/DSCHit.h>


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
jerror_t DMCTrigger_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	// Get attenuation parameters
	double L_over_2 = DBCALGeometry::BCALFIBERLENGTH/2.0;
	double Xo = DBCALGeometry::ATTEN_LENGTH;
	unattenuate_to_center = exp(+L_over_2/Xo);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCTrigger_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// See comments in DMCTrigger_factory.h

	vector<const DBCALHit*> bcalhits;
	vector<const DFCALHit*> fcalhits;
	vector<const DSCHit*> schits;
	loop->Get(bcalhits);
	loop->Get(fcalhits);
	loop->Get(schits);

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
	
	double BCAL_Eupstream = 0.0;
	double BCAL_Edownstream = 0.0;
	for(unsigned int i=0; i< bcalhits.size(); i++){
		const DBCALHit* bcalhit = bcalhits[i];
		
		if(bcalhit->end == DBCALGeometry::kUpstream){
			BCAL_Eupstream += bcalhit->E;
		}else{
			BCAL_Edownstream += bcalhit->E;
		}
	}
	
	// Calculate "energy" sum for BCAL in GeV-ish units
	double Ebcal = unattenuate_to_center * (BCAL_Eupstream + BCAL_Edownstream)/2.0;
	
	// Sum up FCAL energy
	double Efcal = 0.0;
	for(unsigned int i=0; i< fcalhits.size(); i++){
		const DFCALHit* fcalhit = fcalhits[i];
		
		Efcal += fcalhit->E;
	}
	
	// Number of start counter hits
	unsigned int Nschits = schits.size();

	DMCTrigger *trig = new DMCTrigger;

	// BCAL and FCAL
	bool sum_cut = (Ebcal + 4.0*Efcal)>=2.0;
	trig->L1a_fired = sum_cut && Ebcal>0.200 && Efcal>0.030;
	trig->L1b_fired = sum_cut && Ebcal>0.030 && Efcal>0.030 && Nschits>0;
	trig->L1c_fired = Efcal>=0.250;

	trig->Ebcal = Ebcal;
	trig->Efcal = Efcal;
	trig->Nschits = Nschits;
	
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

