#include "JEventProcessor_FCALpedestals.h"
#include <JANA/JApplication.h>
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALHit.h"
#include "FCAL/DFCALDigiHit.h"
#include "FCAL/DFCALCluster.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DVertex.h"
#include "DVector3.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include <TTree.h>
#include "DVector3.h"
#include "PID/DParticleID.h"
#include "GlueX.h"
#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_FCALpedestals());
  }
} // "C"


//------------------
// JEventProcessor_FCALpedestals (Constructor)
//------------------
JEventProcessor_FCALpedestals::JEventProcessor_FCALpedestals()
{

}

//------------------
// ~JEventProcessor_FCALpedestals (Destructor)
//------------------
JEventProcessor_FCALpedestals::~JEventProcessor_FCALpedestals()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FCALpedestals::init(void)
{
  // This is called once at program startup. If you are creating
  // and filling historgrams in this plugin, you should lock the
  // ROOT mutex like this:
  //
  japp->RootWriteLock();

  m_tree = new TTree( "FCALpedestals", "FCAL track pedestals" );

 
  m_tree->Branch( "r", &m_r, "p/I" );
  m_tree->Branch( "c", &m_c, "p/I" );
  m_tree->Branch( "chan", &m_chan, "p/I" );
  m_tree->Branch( "pedestal", &m_pedestal, "p/F" );
 
  
 
  japp->RootUnLock();

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FCALpedestals::brun(JEventLoop *eventLoop, 
					     int32_t runnumber)
{

  // get the FCAL z position from the global geometry interface
  DApplication *dapp = 
    dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  const DGeometry *geom = dapp->GetDGeometry(runnumber);
  if( geom ) {

    geom->GetFCALZ( m_FCALfront );
  }
  else{
      
    cerr << "No geometry accessbile." << endl;
    return RESOURCE_UNAVAILABLE;
  }


  // we need an FCAL Geometry object
  vector< const DFCALGeometry* > geomVec;
  eventLoop->Get( geomVec );

  if( geomVec.size() != 1 ){

    cerr << "No geometry accessbile." << endl;
    return RESOURCE_UNAVAILABLE;
  }

  m_fcalGeom = geomVec[0];

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FCALpedestals::evnt(JEventLoop *eventLoop, 
					     uint64_t eventnumber)
{
  
 
  vector< const DFCALDigiHit*  > digiHits;
   eventLoop->Get( digiHits );
   
  japp->RootWriteLock();
  
   for( vector< const DFCALDigiHit* >::const_iterator dHitItr = digiHits.begin();
       dHitItr != digiHits.end(); ++dHitItr ){
const DFCALDigiHit& dHit = (**dHitItr);
 m_chan = m_fcalGeom->channel( dHit.row, dHit.column );
 m_r = m_fcalGeom->row( dHit.row);  
 m_c = m_fcalGeom->column(dHit.column );   
    
     m_pedestal = dHit.pedestal;
     
    m_tree->Fill();
}

japp->RootUnLock(); 

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FCALpedestals::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FCALpedestals::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}


