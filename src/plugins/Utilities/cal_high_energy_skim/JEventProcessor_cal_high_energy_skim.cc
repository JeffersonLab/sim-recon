// $Id$
//
//    File: JEventProcessor_cal_high_energy_skim.cc
// Created: Mon Dec  1 14:57:11 EST 2014 (copied structure from pi0fcalskim plugin)
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include <math.h>
#include <TLorentzVector.h>
#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include "JEventProcessor_cal_high_energy_skim.h"


#include "TRACKING/DMCThrown.h"
// Routine used to create our JEventProcessor
#include "PID/DVertex.h"
#include "DANA/DApplication.h"
#include "JANA/JApplication.h"
#include "JANA/JFactory.h"
#include "BCAL/DBCALShower.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALCluster.h"
#include "RF/DRFTime.h"
#include "PID/DEventRFBunch.h"
#include "TRIGGER/DTrigger.h"

#include "DLorentzVector.h"
#include "TTree.h"
#include "units.h"
#include "ANALYSIS/DAnalysisUtilities.h"

extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_cal_high_energy_skim());
  }
} // "C"


//------------------
// JEventProcessor_cal_high_energy_skim (Constructor)
//------------------
JEventProcessor_cal_high_energy_skim::JEventProcessor_cal_high_energy_skim()
{

  MIN_BCAL_E = 2.;
  MIN_FCAL_E = 5.;

  //WRITE_EVIO = 1;
  MAKE_DIAGNOSTICS = 0;   // only save ROOT histograms if we need to

  //gPARMS->SetDefaultParameter("CALHIGHENERGY:WRITE_EVIO", WRITE_EVIO );
  gPARMS->SetDefaultParameter("CALHIGHENERGY:DIAGNOSTICS", MAKE_DIAGNOSTICS );
  gPARMS->SetDefaultParameter("CALHIGHENERGY:MIN_BCAL_E" , MIN_BCAL_E );
  gPARMS->SetDefaultParameter("CALHIGHENERGY:MIN_FCAL_E" , MIN_FCAL_E );
  
  num_epics_events = 0;
  
}

//------------------
// ~JEventProcessor_cal_high_energy_skim (Destructor)
//------------------
JEventProcessor_cal_high_energy_skim::~JEventProcessor_cal_high_energy_skim()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_cal_high_energy_skim::init(void)
{
	h_FCAL_shen = new TH1F("FCAL_cluster_energy", "FCAL Cluster Energy; E(cluster), Gev; Counts", 120, 0, 12.);
	h_BCAL_shen = new TH1F("BCAL_shower_energy", "BCAL Shower Energy; E(shower), Gev; Counts", 120, 0, 6.);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_cal_high_energy_skim::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_cal_high_energy_skim::evnt(JEventLoop *loop, uint64_t eventnumber)
{
  vector< const DBCALShower* > locBCALShowers;
  vector< const DFCALShower* > locFCALShowers;
  vector< const DFCALCluster* > locFCALClusters;
  //vector< const DTrackTimeBased*> locTrackTimeBased;
  //loop->Get(locTrackTimeBased);
  //vector<const DVertex*> kinfitVertex;
  //loop->Get(kinfitVertex);

  // initialization
  bool to_save_event = false;
  const DEventWriterEVIO* locEventWriterEVIO = NULL;
  loop->GetSingle(locEventWriterEVIO);
  
  // always write out BOR events
  if(loop->GetJEvent().GetStatusBit(kSTATUS_BOR_EVENT)) {
	  //jout << "Found BOR!" << endl;
	  locEventWriterEVIO->Write_EVIOEvent( loop, "cal_high_energy_skim" );
	  return NOERROR;
  }

  // write out the first few EPICS events to save run number & other meta info
  if(loop->GetJEvent().GetStatusBit(kSTATUS_EPICS_EVENT) && (num_epics_events<5)) {
	  //jout << "Found EPICS!" << endl;
	  locEventWriterEVIO->Write_EVIOEvent( loop, "cal_high_energy_skim" );
	  num_epics_events++;
	  return NOERROR;
  }
  
  // only look at potential physics events
  const DTrigger* locTrigger = NULL;
  loop->GetSingle(locTrigger);
  if(!locTrigger->Get_IsPhysicsEvent())
	  return NOERROR;
  
  // Check to see if we have a high energy BCAL shower
  loop->Get(locBCALShowers);

  if(MAKE_DIAGNOSTICS) japp->RootWriteLock();
  for(unsigned int j = 0; j < locBCALShowers.size(); ++j) {
	  
	  if(MAKE_DIAGNOSTICS) h_BCAL_shen->Fill( locBCALShowers[j]->E_raw );

	  if( locBCALShowers[j]->E_raw > MIN_BCAL_E )
		  to_save_event = true;
	  
	  if(to_save_event)
		  break;
  }
  if(MAKE_DIAGNOSTICS) japp->RootUnLock();


  // if we aren't already going to save these events, save ones with a high energy FCAL shower
  if(!to_save_event || MAKE_DIAGNOSTICS) {

	  loop->Get(locFCALClusters);
	  if(MAKE_DIAGNOSTICS) japp->RootWriteLock();
	  for(unsigned int j = 0; j < locFCALClusters.size(); ++j) {

		  if(MAKE_DIAGNOSTICS) h_FCAL_shen->Fill( locFCALClusters[j]->getEnergy() );
	  
		  if( locFCALClusters[j]->getEnergy() > MIN_FCAL_E )
			  to_save_event = true;
	  
		  if(to_save_event)
			  break;
	  }
	  if(MAKE_DIAGNOSTICS) japp->RootUnLock();
  }

  if(to_save_event) {
	  locEventWriterEVIO->Write_EVIOEvent( loop, "cal_high_energy_skim" );
  }
  
  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_cal_high_energy_skim::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// Fin
//------------------
jerror_t JEventProcessor_cal_high_energy_skim::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}

