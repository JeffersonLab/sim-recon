// $Id$
//
//    File: JEventProcessor_EventTagPi0.cc
// Created: Fri Feb  5 23:23:22 EST 2016
// Creator: davidl (on Darwin harriet 13.4.0 i386)
//

#include "JEventProcessor_EventTagPi0.h"
using namespace jana;

#include <DANA/DStatusBits.h>
#include <FCAL/DFCALHit.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_EventTagPi0());
}
} // "C"


//------------------
// JEventProcessor_EventTagPi0 (Constructor)
//------------------
JEventProcessor_EventTagPi0::JEventProcessor_EventTagPi0()
{

}

//------------------
// ~JEventProcessor_EventTagPi0 (Destructor)
//------------------
JEventProcessor_EventTagPi0::~JEventProcessor_EventTagPi0()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_EventTagPi0::init(void)
{
	Emin_MeV = 200.0;
	Rmin_cm  = 30.0;
	
	gPARMS->SetDefaultParameter("PI0TAG:Emin_MeV", Emin_MeV , "Minimum energy in MeV of each single block hit to tag event as FCAL pi0");
	gPARMS->SetDefaultParameter("PI0TAG:Rmin_cm" , Rmin_cm  , "Minimum distance in cm between single blocks with energy > PI0TAG:Emin_MeV to tag event as FCAL pi0");

	Rmin_cm_2 = Rmin_cm*Rmin_cm;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_EventTagPi0::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_EventTagPi0::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	vector<const DFCALHit*> fcalhits;
	loop->Get(fcalhits);
	
	if(fcalhits.size() < 2) return NOERROR;
	
	for(uint32_t i=0; i<(fcalhits.size()-1); i++){
		const DFCALHit *hit1 = fcalhits[i];
		if( hit1->E < Emin_MeV ) continue;

		for(uint32_t j=i+1; j<fcalhits.size(); j++){
			const DFCALHit *hit2 = fcalhits[j];
			if( hit2->E < Emin_MeV ) continue;
	
			double deltaX = hit1->x - hit2->x;
			double deltaY = hit1->y - hit2->y;
			double r2 = deltaX*deltaX + deltaY*deltaY;
			if( r2 <= Rmin_cm_2 ){
				JEvent &jevent = loop->GetJEvent();
				jevent.SetStatusBit(kSTATUS_FCAL_PI0);
				jevent.SetStatusBit(kSTATUS_PI0);
				return NOERROR;
			}
		}
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_EventTagPi0::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_EventTagPi0::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

