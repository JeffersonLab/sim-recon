// $Id$
//
//    File: JEventProcessor_BCAL_SiPM_saturation.cc
//          Modified file from BCAL_neutron_discriminator.cc   ES  5/10/2018
// Created: Thu Apr  5 16:36:00 EDT 2018
// Creator: dalton (on Linux gluon119.jlab.org 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_BCAL_SiPM_saturation.h"
#include "HistogramTools.h"
#include "TRACKING/DMCThrown.h"
// #include "BCAL/DBCALHit.h"
// #include "BCAL/DBCALTDCHit.h"
// #include "BCAL/DBCALCluster.h"
// #include "BCAL/DBCALDigiHit.h"
#include "BCAL/DBCALPoint.h"
// #include "BCAL/DBCALUnifiedHit.h"
// #include "BCAL/DBCALGeometry.h"
#include "BCAL/DBCALShower.h"
#include "DANA/DStatusBits.h"
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

    attenuation_parameters.clear();
    eventLoop->GetCalib("/BCAL/attenuation_parameters", attenuation_parameters);
    /*int channel = 0;
       for (int module=1; module<=48; module++) {
           for (int layer=1; layer<=4; layer++) {
               for (int sector=1; sector<=4; sector++) {
                   printf("%2i  %2i  %2i %12.4f %12.4f %12.4f\n",
                          module,layer,sector,
                          attenuation_parameters[channel][0],
                          attenuation_parameters[channel][1],
                          attenuation_parameters[channel][2]);
               }
               channel++;
           }
	   }*/

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_SiPM_saturation::evnt(JEventLoop *loop, uint64_t eventnumber)
{
 

  // Check to see if this is a physics event
  if (loop->GetJEvent().GetStatusBit(kSTATUS_PHYSICS_EVENT)) {

  // Check if thrown information is available
  vector<const DMCThrown*> MCThrowns;
  loop->Get(MCThrowns);
  unsigned int NumThrown=MCThrowns.size();
  double Ethrown=0;
  if (NumThrown == 1) {
    const DMCThrown* locMCThrown = MCThrowns[0];
    Ethrown = locMCThrown->energy();
  }

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
        float Eshower = locBCALShower->E;
        /*float E_preshower = locBCALShower->E_preshower;
        float z = locBCALShower->z;
        float x = locBCALShower->x;
        float y = locBCALShower->y;
        float R = sqrt(x*x+y*y);
        float sigLong = locBCALShower->sigLong;    // rho
        float sigTrans = locBCALShower->sigTrans;  // phi
        float sigTheta = locBCALShower->sigTheta;*/

        // cout << " Shower i=" << i << " E_shower=" << Eshower << " E_preshower=" << E_preshower << " x=" << x << " y=" << y << " z=" << z << " R=" << R << " sigLong=" << sigLong << " sigTrans=" << sigTrans << " sigTheta=" << sigTheta << endl;

	Int_t nbins=100;

	// Fill histogram for showers
        Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Ethrown", Ethrown,
                         "BCAL SiPM Saturation; Thrown Energy (GeV)",4*nbins,0,10);

        Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Eshower", Eshower,
                         "BCAL SiPM Saturation; Shower Energy (GeV)",4*nbins,0,10);

        // Fill 2D histograms
        Fill2DHistogram ("BCAL_SiPM_saturation", "Hists2D", "Eshower_vs_Ethrown", Ethrown, Eshower,
                         "BCAL SiPM Saturation; Shower Energy (GeV); Thrown Energy (GeV)",
                         4*nbins,0,10,4*nbins,0,10);
        Fill2DHistogram ("BCAL_SiPM_saturation", "Hists2D", "EDiff_vs_Ethrown", Ethrown, Eshower - Ethrown,
                         "BCAL SiPM Saturation; Thrown Energy (GeV); (Shower - Thrown) Energy (GeV)",
                         4*nbins,0,10,nbins,-0.5,0.5);
        Fill2DHistogram ("BCAL_SiPM_saturation", "Hists2D", "EDiff/Ethrown_vs_Ethrown", Ethrown, (Eshower - Ethrown)/Ethrown,
                         "BCAL SiPM Saturation; Thrown Energy (GeV); (EShower - EThrown)/Ethrown",
                         4*nbins,0,10,nbins,-0.1,0.1);
	
        // Get vector of points in this shower
        vector<const DBCALPoint*> Points;
		locNeutralShower->Get(Points);
        uint Ncell = Points.size();

        for (unsigned int j = 0; j < Ncell; j++){
            const DBCALPoint* locPoint = Points[j];
            float t = locPoint->t();
            float z = locPoint->z();
	    float Ept = locPoint->E();
	    int module= locPoint->module();
	    int layer = locPoint->layer();
	    int sector = locPoint->sector();
	    // cout << " j=" << j << " t=" << t << endl;

	    // calculate point energy (code taken from DBCALPoint.cc)
	    float z_bcal_center = 212;
	    float fibLen = 390;
	    float zLocal = z - z_bcal_center; 


	    float dUp = 0.5 * fibLen + zLocal;
	    float dDown = 0.5 * fibLen - zLocal;
	    if (dUp>fibLen)   dUp=fibLen;
	    if (dUp<0)        dUp=0;
	    if (dDown>fibLen) dDown=fibLen;
	    if (dDown<0)      dDown=0;
	    int channel = (module-1)*16 + (layer-1)*4 + (sector-1);
	    // cout << " module=" << module << " layer=" << layer << " sector=" << sector << " channel=" << channel << endl;
	    // cout << " attenuation length 00=" << attenuation_parameters[channel][0] << endl;
	    float attenuation_length = attenuation_parameters[channel][0];
	    float attUp = exp( -dUp / attenuation_length );
	    float attDown = exp( -dDown / attenuation_length );
 

	    if (VERBOSE>=3) cout << " VERBOSE >=3" << " t=" << t << " z=" << z << endl;


        // Fill 1D histograms

        Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "NCell", Ncell,
                         "BCAL SiPM Saturation; Number of cells",nbins,0,100);

        Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "layer", layer,
                         "BCAL SiPM Saturation; Layer Number",5,0,5);

        Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Ept", Ept,
                         "BCAL SiPM Saturation; Point Energy (GeV)",4*nbins,0,10);

	// cout << " Point: Ept=" << Ept << endl;
	float upHit=0;
	float downHit=0;
	float uppeak=0;
	float downpeak=0;

             vector<const DBCALHit*> Hits;
		locPoint->Get(Hits);
                uint Nhits = Hits.size();

             for (unsigned int j = 0; j < Nhits; j++){
                const DBCALHit* locHit = Hits[j];
	        float Ehit = locHit->E;
	        /*int module = locHit->module;
	        int layer = locHit->layer;
	        int sector = locHit->sector;*/
	        int end = locHit->end;
	        int pulse_peak = locHit->pulse_peak;

		// cout << " module=" << module << " layer=" << layer << " sector=" << sector << " pulse_peak/Ehit=" << pulse_peak/Ehit << endl;

		if (end == 0) {
		  upHit = Ehit;
		  uppeak = pulse_peak;
		}
		if (end == 1) {
		  downHit = Ehit;
		  downpeak = pulse_peak;
		}


	    nbins=5100;

            // Fill 1D histograms
	    if (layer == 1) {
            Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Hit pulse_peak layer=1", pulse_peak,
                         "BCAL SiPM Saturation; Hit Pulse_peak (counts)",nbins,-100,5000);

            Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Hit integral layer=1", Ehit,
			     "BCAL SiPM Saturation; Hit integral (counts)",nbins,0, 10);
	    }
	    else if (layer == 2) {
            Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Hit pulse_peak layer=2", pulse_peak,
                         "BCAL SiPM Saturation; Hit Pulse_peak (counts)",nbins,-100,5000);

            Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Hit integral layer=2", Ehit,
			     "BCAL SiPM Saturation; Hit integral (counts)",nbins,0, 10);
	    }
	    else if (layer == 3) {
            Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Hit pulse_peak layer=3", pulse_peak,
                         "BCAL SiPM Saturation; Hit Pulse_peak (counts)",nbins,-100,5000);

            Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Hit integral layer=3", Ehit,
			     "BCAL SiPM Saturation; Hit integral (counts)",nbins,0, 10);
	    }
	    else if (layer == 4) {
            Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Hit pulse_peak layer=4", pulse_peak,
                         "BCAL SiPM Saturation; Hit Pulse_peak (counts)",nbins,-100,5000);

            Fill1DHistogram ("BCAL_SiPM_saturation", "Hists1D", "Hit integral layer=4", Ehit,
			     "BCAL SiPM Saturation; Hit integral (counts)",nbins,0, 10);
	    }
	    else {
	      cout << " ***Illegal layer=" << layer << endl;
	    }


	     }



		// use these to correct the energy
		float E_US =  ( upHit / attUp );
		float E_DS =  ( downHit / attDown );
		float Ecalc =  ( E_US + E_DS ) / 2;

		/*cout << " i=" << i << " module=" << module <<  " layer=" <<  layer << " sector=" << sector 
                     << " lambda=" << attenuation_length << " Eup=" << upHit << " Edown=" << downHit << " Point: Ept=" 
                     << Ept << " Ecalc=" << Ecalc << " Diff=" << Ept-Ecalc <<  endl;*/

	}



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

