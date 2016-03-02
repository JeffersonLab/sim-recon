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

// Set number of channels
const int nChan = 2800;

// Define Histograms
static TH1F* pedestal[nChan];
//static TH1F* timingoffset_vert[nChan];



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

  //TDirectory *main1 = gDirectory;
  //gDirectory->mkdir("horizontal_offsets")->cd();
  for (int i = 0; i < nChan; ++i) {
    pedestal[i] = new TH1F(Form("pedestal_%i",i),Form("Pedestal for Channel %i",i),500,90,120);
  }
 
  
 
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
 //m_r = m_fcalGeom->row( dHit.row);  
 //m_c = m_fcalGeom->column(dHit.column );   
    
     m_pedestal = dHit.pedestal;
     
     //cout << "Channel: " << m_chan << " Pedestal: " << m_pedestal << endl;
     pedestal[m_chan]->Fill(m_pedestal);  
   // m_tree->Fill();
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


