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
#include <map>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <TDirectory.h>
#include <TH1I.h>

using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

//FCAL has only 2800 active channels but and a bunch of inactive blocks. Redundancy in number of channels is to make sure there are enough histograms available to fill all active channels.
const int nChan = 5000;

// Define Histograms
static TH1I* pedestal[nChan][nChan];
//static TH1I* pedestal[nChan];
//static TH1F* timingoffset_vert[nChan];
//  map < int, pair <int,int> > mapchan;


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

  TDirectory *main = gDirectory;
  gDirectory->mkdir("FCAL_pedestals")->cd();
  
  /*
 for (int i = 0; i < nChan; ++i) {
    pedestal[i] = new TH1I(Form("pedestal_%i",i),Form("Pedestal for Channel %i",i),500,90,120);
  }
 */
  for (int i = 0; i <59; ++i) {
  	for (int j = 0; j <59; ++j) {
   		 pedestal[i][j] = new TH1I(Form("Pedestal_%i_%i",i,j),Form("Pedestal for Channel: %i_%i",i,j),500,90,120);
 		}
	}
 main->cd();
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

 m_r = m_fcalGeom->row( dHit.row) ;  
 m_c = m_fcalGeom->column(dHit.column ) ;   
     /*
       if( !m_fcalGeom->isBlockActive( m_r, m_c ) )
       {
       cout << "Invalid Block (X,Y): " << "\t" << m_c-29 << "\t" << m_r-29 << "\t" << m_fcalGeom->channel( dHit.row, dHit.column ) << endl;
       }
       */
              if( !m_fcalGeom->isBlockActive( m_r, m_c ) ) continue;
       //  m_chan = m_fcalGeom->channel( dHit.row, dHit.column );
     m_pedestal = dHit.pedestal;
            if( m_pedestal > 0 ) {
     //cout << "Channel: " << m_chan << " Pedestal: " << m_pedestal << endl;
 // cout << "Channel: " << m_chan << " Pedestal: " << m_pedestal << " X: "<< m_c - 29 << " Y: " << m_r - 29 << endl;   
 //pedestal[m_chan] = new TH1F(Form("Pedestal_%i_%i",m_r,m_c),Form("Pedestal for Channel: %i_%i",m_r,m_c),500,90,120);
// if(m_pedestal > 0) {
     pedestal[m_r][m_c]->Fill(m_pedestal);  
  //   }
      //pedestal[m_chan]->Fill(m_pedestal); 
      //SetName(Form("Pedestal_%i_%i",m_c,m_r));
   // m_tree->Fill();
   
//   mapchan[m_chan] = std::make_pair(m_c,m_r);
   }
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

  // for (int i = 0; i < mapchan.size(); ++i) {
//cout << mapchan->first[i] << "\t" << mapchan->second[i] << endl;
//  }
 /*
 for (map < int, pair <int,int> >::const_iterator it = mapchan.begin(); it != mapchan.end(); it++) {
    std::cout << "Key: " << it->first << "\t"<< "Values: " << it->second.first << "\t" << it->second.second << endl;
}

   cout << "Size of map: " << mapchan.size() << endl;
   
   */
  return NOERROR;
}


