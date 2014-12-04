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

  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MIN_MASS", MIN_MASS );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MAX_MASS", MAX_MASS );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MIN_E", MIN_E );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MIN_R", MIN_R );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MAX_DT", MAX_DT );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MAX_ETOT", MAX_ETOT );
  gPARMS->SetDefaultParameter( "PI0FCALSKIM:MIN_BLOCKS", MIN_BLOCKS );
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

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_pi0fcalskim::brun(JEventLoop *eventLoop, int runnumber)
{

  eventLoop->GetSingle(dEventWriterEVIO);

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_pi0fcalskim::evnt(JEventLoop *loop, int eventnumber)
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

    dEventWriterEVIO->Write_EVIOEvent( loop, "pi0fcalskim" );
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
// fini
//------------------
jerror_t JEventProcessor_pi0fcalskim::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}

