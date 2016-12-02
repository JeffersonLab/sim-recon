// $Id$
//
//    File: JEventProcessor_FCALpulsepeak.cc
// Created: Tue Sep 27 11:18:28 EDT 2016
// Creator: asubedi (on Linux stanley.physics.indiana.edu 2.6.32-573.18.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_FCALpulsepeak.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <stdint.h>
#include <vector>
#include <iostream>
#include "FCAL/DFCALHit.h"
#include "FCAL/DFCALDigiHit.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALCluster.h"
#include "FCAL/DFCALShower.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250PulsePedestal.h"
#include "DAQ/Df250PulseData.h"
#include "units.h"
#include "DLorentzVector.h"
#include "DVector3.h"
#include "HDGEOMETRY/DGeometry.h"
#include "DANA/DApplication.h"
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <TDirectory.h>
#include <TH1I.h>
#include <TH2F.h>
#include "HistogramTools.h"


const int nChan = 2800;

// Define Histograms
static TH1I* pulsepeak[nChan];


extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_FCALpulsepeak());
}
} // "C"


//------------------
// JEventProcessor_FCALpulsepeak (Constructor)
//------------------
JEventProcessor_FCALpulsepeak::JEventProcessor_FCALpulsepeak()
{

}

//------------------
// ~JEventProcessor_FCALpulsepeak (Destructor)
//------------------
JEventProcessor_FCALpulsepeak::~JEventProcessor_FCALpulsepeak()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FCALpulsepeak::init(void)
{

	//This is called once at program startup. 
TDirectory *main = gDirectory;
  gDirectory->mkdir("FCAL_pulsepeak")->cd();
  
  
  for (int i = 0; i < nChan; ++i) {
    pulsepeak[i] = new TH1I(Form("peak_%i",i),Form("Pulsepeak for Channel %i",i),500,10,-10);
  }

  main->cd();

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FCALpulsepeak::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}
//------------------
// evnt
//------------------
jerror_t JEventProcessor_FCALpulsepeak::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{

vector< const DFCALDigiHit*  > digiHits;
   eventLoop->Get( digiHits );
   vector< const DFCALGeometry* > geomVec;
    eventLoop->Get( geomVec );
  
  const DFCALGeometry& fcalGeom = *(geomVec[0]);
   
for( vector< const DFCALDigiHit* >::const_iterator dHitItr = digiHits.begin();
       dHitItr != digiHits.end(); ++dHitItr ){
const DFCALDigiHit& dHit = (**dHitItr);
	int m_chan = fcalGeom.channel( dHit.row, dHit.column );
     //uint32_t m_time = dHit.pulse_time;
     uint32_t m_peak = dHit.pulse_peak;
     uint32_t m_pedestal = dHit.pedestal;
     //uint32_t m_integral = dHit.pulse_integral;
     
      japp->RootFillLock(this);
     if (m_peak > 0 && m_pedestal > 95){
     
      pulsepeak[m_chan]->Fill(m_peak); 
      
      }

	japp->RootFillUnLock(this);
}


/*
	  vector< const DFCALGeometry* > geomVec;
  vector< const DFCALDigiHit*  > digiHits;
  vector< const DFCALHit*      > hits;
    eventLoop->Get( geomVec );
  eventLoop->Get( digiHits );
  eventLoop->Get( hits );
  
  
  const DFCALGeometry& fcalGeom = *(geomVec[0]);
 

map< const DFCALDigiHit*, const Df250PulseData* > pd_cache;

  for( vector< const DFCALDigiHit* >::const_iterator dHitItr = digiHits.begin();
       dHitItr != digiHits.end(); ++dHitItr ){
  
    // fetch lower level FADC data
    const Df250PulseData*     pulseDat = NULL;

    const DFCALDigiHit& dHit = (**dHitItr);
    dHit.GetSingle( pulseDat );
    pd_cache[&dHit] = pulseDat;
  }
  
  for( vector< const DFCALDigiHit* >::const_iterator dHitItr = digiHits.begin();
       dHitItr != digiHits.end(); ++dHitItr ){

    const DFCALDigiHit& dHit = (**dHitItr);
   
    const Df250PulseData* pulseDat = pd_cache[&dHit];
    
    
    if( pulseDat){
      
      // m_x = dHit.column  - 29;
      // m_y = dHit.row  - 29;
       m_chan = fcalGeom.channel( dHit.row, dHit.column );
    //cout << "pedestal: " << pulseDat->pedestal << " peak: " << pulseDat->pulse_peak << endl;
    
       m_peak = pulseDat->pulse_peak;
   
       japp->RootFillLock(this);
       if (m_peak > 0){
       // cout << "Channel: " << m_chan << " peak: " << m_peak << endl;
      pulsepeak[m_chan]->Fill(m_peak); 
   // Fill the only important histogram
                //Fill2DHistogram("FCALpulsepeak", "", "FCAL Pulse Peak Vs. Channel Number",
                        //m_chan, m_peak,
                        //"Pulse Peak ;Channel Number; Pulse Peak",
                        //2800,-0.5,2799.5,4000,100, 4100); // Channels are numbered from zero...AGAINST CONVENTION    
     }
japp->RootFillUnLock(this);
 

	}
     }
*/  
  
  





	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FCALpulsepeak::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FCALpulsepeak::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

