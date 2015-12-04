// $Id$
//
//    File: JEventProcessor_pi0fcalskim.cc
// Created: Mon Dec  1 14:57:11 EST 2014
// Creator: shepherd (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include <math.h>

#include "JEventProcessor_pi0fcalskim.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include "JANA/JApplication.h"
#include "JANA/JFactory.h"
#include "FCAL/DFCALCluster.h"
#include "DLorentzVector.h"
#include "TTree.h"
#include "units.h"

extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_pi0fcalskim());
  }
} // "C"


//------------------
// JEventProcessor_pi0fcalskim (Constructor)
//------------------
JEventProcessor_pi0fcalskim::JEventProcessor_pi0fcalskim()
{

  MIN_MASS   = 0.03; // GeV
  MAX_MASS   = 0.30; // GeV
  MIN_E      =  1.0; // GeV (photon energy cut)
  MIN_R      =   20; // cm  (cluster distance to beam line)
  MAX_DT     =   10; // ns  (cluster time diff. cut)
  MAX_ETOT   =   12; // GeV (max total FCAL energy)
  MIN_BLOCKS =    2; // minumum blocks per cluster

  WRITE_ROOT = 0;
  WRITE_EVIO = 1;

  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MIN_MASS", MIN_MASS );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MAX_MASS", MAX_MASS );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MIN_E", MIN_E );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MIN_R", MIN_R );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MAX_DT", MAX_DT );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MAX_ETOT", MAX_ETOT );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MIN_BLOCKS", MIN_BLOCKS );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:WRITE_ROOT", WRITE_ROOT );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:WRITE_EVIO", WRITE_EVIO );
}

//------------------
// ~JEventProcessor_pi0fcalskim (Destructor)
//------------------
JEventProcessor_pi0fcalskim::~JEventProcessor_pi0fcalskim()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_pi0fcalskim::init(void)
{
  dEventWriterEVIO = NULL;

  if( ! ( WRITE_ROOT || WRITE_EVIO ) ){

    cerr << "No output mechanism has been specified." << endl;
    return UNRECOVERABLE_ERROR;
  }

  if( WRITE_ROOT ){

    japp->RootWriteLock();

    m_tree = new TTree( "cluster", "Cluster Tree for Pi0 Calibration" );
    m_tree->Branch( "nClus", &m_nClus, "nClus/I" );
    m_tree->Branch( "hit0", m_hit0, "hit0[nClus]/I" );
    m_tree->Branch( "px", m_px, "px[nClus]/F" );
    m_tree->Branch( "py", m_py, "py[nClus]/F" );
    m_tree->Branch( "pz", m_pz, "pz[nClus]/F" );

    m_tree->Branch( "nHit", &m_nHit, "nHit/I" );
    m_tree->Branch( "chan", m_chan, "chan[nHit]/I" );
    m_tree->Branch( "e", m_e, "e[nHit]/F" );

    japp->RootUnLock();
  }

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_pi0fcalskim::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  eventLoop->GetSingle(dEventWriterEVIO);

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_pi0fcalskim::evnt(JEventLoop *loop, uint64_t eventnumber)
{
  vector< const DFCALCluster* > clusterVec;
  loop->Get( clusterVec );

  if( clusterVec.size() < 2 ) return NOERROR;

  bool hasCandidate = false;
  double eTot = 0;

  for( vector< const DFCALCluster*>::const_iterator clus1Itr = clusterVec.begin();
       clus1Itr != clusterVec.end(); ++clus1Itr ){

    eTot += (**clus1Itr).getEnergy();

    for( vector< const DFCALCluster*>::const_iterator clus2Itr = clus1Itr + 1;
	 clus2Itr != clusterVec.end(); ++clus2Itr ){

      const DFCALCluster& clusL = 
	( (**clus1Itr).getEnergy() > (**clus2Itr).getEnergy() ? 
	  (**clus2Itr) : (**clus1Itr) );

      const DFCALCluster& clusH = 
	( (**clus1Itr).getEnergy() > (**clus2Itr).getEnergy() ? 
	  (**clus1Itr) : (**clus2Itr) );

      double clusLX = clusL.getCentroid().X();
      double clusLY = clusL.getCentroid().Y();
      double rL = sqrt( clusLX * clusLX + clusLY * clusLY );
      double eL = clusL.getEnergy();
      double tL = clusL.getTime();
      int nHitL = clusL.GetHits().size();

      double clusHX = clusH.getCentroid().X();
      double clusHY = clusH.getCentroid().Y();
      double rH = sqrt( clusHX * clusHX + clusHY * clusHY );
      double eH = clusH.getEnergy();
      double tH = clusH.getTime();
      int nHitH = clusH.GetHits().size();

      DVector3 clusLMom = clusL.getCentroid(); 
      clusLMom.SetMag( eL );

      DVector3 clusHMom = clusH.getCentroid(); 
      clusHMom.SetMag( eH );
    
      double dt = fabs( tL - tH );

      DLorentzVector gamL( clusLMom, clusLMom.Mag() );
      DLorentzVector gamH( clusHMom, clusHMom.Mag() );

      double mass = ( gamL + gamH ).M();

      hasCandidate |= 
	( ( eL > MIN_E ) &&
	  ( dt < MAX_DT ) &&
	  ( rL > MIN_R ) && ( rH > MIN_R ) &&
	  ( nHitL >= MIN_BLOCKS ) && ( nHitH >= MIN_BLOCKS ) &&
	  ( mass > MIN_MASS ) && ( mass < MAX_MASS  ) );
    }
  }

  if( hasCandidate && ( eTot < MAX_ETOT ) ){

    if( WRITE_EVIO ){

      dEventWriterEVIO->Write_EVIOEvent( loop, "pi0fcalskim" );
    }

    if( WRITE_ROOT ){

      japp->RootWriteLock();
      writeClustersToRoot( clusterVec );
      japp->RootUnLock();
    }
  }

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_pi0fcalskim::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// Fin
//------------------
jerror_t JEventProcessor_pi0fcalskim::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}

void 
JEventProcessor_pi0fcalskim::writeClustersToRoot( const vector< const DFCALCluster* > clusVec ){

  // this code must run serially -- obtain a lock before
  // entering this function
 
  m_nHit = 0;
  m_nClus = 0;

  // hit and cluster indices
  int& iH = m_nHit;
  int& iC = m_nClus;
  
  for( vector< const DFCALCluster* >::const_iterator clusItr = clusVec.begin();
       clusItr != clusVec.end(); ++clusItr ){

    // if we exceed max clusters abort writing for this event
    if( iC == kMaxClus ) return;

    const DFCALCluster& clus = (**clusItr);

    if( ( clus.getCentroid().Perp() < 20*k_cm ) ||
	( clus.getEnergy() < 1*k_GeV ) ||
	( clus.GetHits().size() < 2 ) ) continue;

    DVector3 gamMom = clus.getCentroid(); 
    gamMom.SetMag( clus.getEnergy() );

    m_hit0[iC] = iH;
    m_px[iC] = gamMom.X();
    m_py[iC] = gamMom.Y();
    m_pz[iC] = gamMom.Z();
    
    const vector<DFCALCluster::DFCALClusterHit_t>& hits = clus.GetHits();

    for( unsigned int i = 0; i < hits.size(); ++i ){

      // if we exceed max hits abort this event and return
      if( iH == kMaxHits ) return;

      m_chan[iH] = hits[i].ch;
      m_e[iH] = hits[i].E;
      ++iH;
    }
    ++iC;
  }

  m_tree->Fill();
}
