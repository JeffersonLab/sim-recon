// $Id$
//
//    File: JEventProcessor_BCAL_SiPM_saturation.cc
//          Modified file from BCAL_neutron_discriminator.cc   ES  5/10/2018
// Created: Thu Apr  5 16:36:00 EDT 2018
// Creator: dalton (on Linux gluon119.jlab.org 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_BCAL_SiPM_saturation.h"
#include "HistogramTools.h"
// #include "BCAL/DBCALHit.h"
// #include "BCAL/DBCALTDCHit.h"
// #include "BCAL/DBCALCluster.h"
// #include "BCAL/DBCALDigiHit.h"
#include "BCAL/DBCALPoint.h"
// #include "BCAL/DBCALUnifiedHit.h"
// #include "BCAL/DBCALGeometry.h"
#include "BCAL/DBCALShower.h"
// #include "DANA/DStatusBits.h"
// #include "PID/DChargedTrack.h"
// #include "PID/DEventRFBunch.h"
// #include "PID/DDetectorMatches.h"
#include "PID/DNeutralShower.h"
// #include "PID/DVertex.h"
// #include "TRACKING/DTrackTimeBased.h"
// #include "TRIGGER/DL1Trigger.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_BCAL_SiPM_saturation());
}
} // "C"


//------------------
// JEventProcessor_BCAL_SiPM_saturation (Constructor)
//------------------
JEventProcessor_BCAL_SiPM_saturation::JEventProcessor_BCAL_SiPM_saturation()
{
  	VERBOSE = 0;
	if(gPARMS){
		gPARMS->SetDefaultParameter("BCAL_SiPM_saturation:VERBOSE", VERBOSE, "Verbosity level");
	}
}

//------------------
// ~JEventProcessor_BCAL_SiPM_saturation (Destructor)
//------------------
JEventProcessor_BCAL_SiPM_saturation::~JEventProcessor_BCAL_SiPM_saturation()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_SiPM_saturation::init(void)
{
	// This is called once at program startup. 

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_SiPM_saturation::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_SiPM_saturation::evnt(JEventLoop *loop, uint64_t eventnumber)
{

    // Get vector of neutral showers in event
  	vector<const DNeutralShower*> NeutralShowers;
  	loop->Get(NeutralShowers);

    unsigned int NumShowers=NeutralShowers.size();

    // loop over neutral showers
    for (unsigned int i = 0; i < NumShowers; i++){
        const DNeutralShower* locNeutralShower = NeutralShowers[i];
        // Must be BCAL shower
        DetectorSystem_t locDetector = locNeutralShower->dDetectorSystem;
        if (locDetector != SYS_BCAL) continue;
        // Get shower properties
        vector<const DBCALShower*> BCALShowers;
	    locNeutralShower->Get(BCALShowers);
        const DBCALShower* locBCALShower = BCALShowers[0];
        float E = locBCALShower->E;
        float E_preshower = locBCALShower->E_preshower;
        float z = locBCALShower->z;
        float x = locBCALShower->x;
        float y = locBCALShower->y;
        float R = sqrt(x*x+y*y);
        float sigLong = locBCALShower->sigLong;    // rho
        float sigTrans = locBCALShower->sigTrans;  // phi
        float sigTheta = locBCALShower->sigTheta;

        cout << " Shower i=" << i << " E=" << E << " E_preshower=" << E_preshower << " x=" << x << " y=" << y << " z=" << z << " R=" << R << " sigLong=" << sigLong << " sigTrans=" << sigTrans << " sigTheta=" << sigTheta << endl;
	
        // Get vector of points in this shower
        vector<const DBCALPoint*> Points;
		locNeutralShower->Get(Points);
        uint Ncell = Points.size();

        for (unsigned int j = 0; j < Ncell; j++){
            const DBCALPoint* locPoint = Points[j];
            float t = locPoint->t();
	    cout << " j=" << j << " t=" << t << endl;

	    if (VERBOSE>=3) cout << " VERBOSE >=3" << " t=" << t << endl;

	Int_t nbins=100;

        // Fill 1D histograms
        Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "NCell", Ncell,
                         "BCAL SiPM Saturation; Number of cells",nbins,0,100);

        // Fill 2D histograms
        Fill2DHistogram ("BCAL_SiPM_saturation", "Hists2D", "RMSt_NCell", Ncell, t,
                         "BCAL SiPM Saturation; Number of cells; RMS of point times  (ns)",
                         nbins,0,100,nbins,0,100);
	}
    }
    return NOERROR;
    }

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_SiPM_saturation::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_SiPM_saturation::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

