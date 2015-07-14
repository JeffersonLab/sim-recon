// $Id$
//
//    File: JEventProcessor_FCAL_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#include <stdint.h>
#include <vector>

#include "JEventProcessor_FCAL_online.h"
#include <JANA/JApplication.h>

#include "FCAL/DFCALHit.h"
#include "FCAL/DFCALDigiHit.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALCluster.h"
#include "FCAL/DFCALShower.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250PulsePedestal.h"
#include "units.h"
#include "DLorentzVector.h"
#include "DVector3.h"
#include "HDGEOMETRY/DGeometry.h"
#include "DANA/DApplication.h"

#include <TDirectory.h>
#include <TH2F.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TProfile.h>

using namespace std;
using namespace jana;


//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_FCAL_online());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_FCAL_online::JEventProcessor_FCAL_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_FCAL_online::~JEventProcessor_FCAL_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_FCAL_online::init(void) {

  // lock all root operations
  japp->RootWriteLock();

  // create root folder for fcal and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("fcal")->cd();

  //
  // these with "dig" use raw channel-level output of the FADC
  //
  fcal_num_events = new TH1D("fcal_num_events","FCAL number of events",1, 0.5, 1.5);

  m_digInt = new TH1I( "digHitE", "FCAL Raw Pulse Integral; Integral [2 V * 4 ns / 4096]; Pulses / ( 100 * 2 V * 4 ns / 4096 )", 300, 0, 30000 );
  m_digCoarseT = new TH1I( "digCoarseTime", "FCAL Raw Pulse Coarse Time; Time [4 ns]; Pulses / 4 ns", 100, 0, 100 );
  m_digCoarseTChan = new TProfile( "digCoarseTChan", "FCAL Coarse Time vs. Channel; Channel Index; Average Time [4 ns]",
				   2800,-0.5, 2799.5 );
  m_digPreciseT = new TH1I( "digPreciseT", "FCAL Raw Pulse Precise Time; Time [62.5 ps]; Pulses / 62.5 ps", 64, -0.5, 63.5 );
  m_digPreciseTChan = new TProfile( "digPreciseTChan", "FCAL Precise Time vs. Channel; Channel Index; Average Time [62.5 ps]",
				   2800,-0.5, 2799.5 );
  m_digT = new TH1I( "digT", "FCAL Pulse Time; t [62.5 ps]; Pulses / 625 ps", 500, 0, 5000 );
  m_digT0 = new TH1I( "digT0", "FCAL Pulse Energy Weighted Average Time; t_{0} [62.5 ps]; Events / 625 ps", 500, 0, 5000 );
  m_digTmT0 = new TH1I( "digTmT0", "FCAL Pulse Local Time; t - t_{0} [62.5 ps]; Pulses / 625 ps", 200, -1000, 1000 );
  m_digTmT02D = new TH2F( "digTmT02D", "FCAL Pulse Local Time [62.5 ps]; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );
  m_digPed = new TH1I( "digPed", "FCAL Pedestal (All Channels); ADC Counts; Pules / 10 Counts", 100, 0, 1000 );
  m_digPedChan = new TProfile( "digPedChan", "FCAL Pedestal vs. Channel; Channel Index; Average Pedestal [ADC Counts]", 2800,-0.5,2799.5 );
  m_digPed2D = new TH2F( "digPed2D", "FCAL Pedestal [ADC Counts]; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );
  m_digPedSq2D = new TH2F( "digPedSq2D", "FCAL Pedestal^{2} [ADC Counts^{2}]; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );
  m_digQual = new TH1I( "digQual", "FCAL Hit Quality; Quality Factor Index; Number of Pulses", 16, -0.5, 15.5 );
  m_digNUnder = new TH1I( "digNUnder", "FCAL Number of Underflows per Event; Number of Underflows; Events / 1", 
			  100, -0.5, 99.5 );
  m_digNOver = new TH1I( "digNOver", "FCAL Number of Overflows per Event; Number of Overflows; Events / 1", 
			 100, -0.5, 99.5 );
  m_digN = new TH1I( "digN", "FCAL Number of Raw Hits per Event; Number of Hits; Events / 5", 600, 0, 3000 );
  m_digPeakV = new TH1I( "digPeakV", "FCAL Pulse Peak; Peak Sample [Volts]; Pulses / 10 mV", 210, 0, 2.1 );
  m_digPeakV2D = new TH2F( "digPeakV2D", "FCAL Pulse Peak [Volts]; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );
  m_digOcc2D = new TH2F( "digOcc2D", "FCAL Pulse Occupancy; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );
  m_digIntVsPeak = new TH2I( "digIntVsPeak", "FCAL Pulse Integral vs. Peak Sample; Peak Sample [ADC Units]; Integral [ADC Units]", 200, 0, 4000, 200, 0, 40000 );
  m_digIntToPeak = new TH1I( "digIntToPeak", "Pedestal Subtracted Integral to Peak Ratio (Peak > 200); Ratio; Pulses / 0.1", 200, 0, 20 );

  //
  // these with "hit" are block level plots after calibration
  //

  m_hitN = new TH1I( "hitN", "FCAL Number of Hits; Number of Hits; Events / 5 Hits", 600, 0, 3000 );
  m_hitE = new TH1I( "hitE", "FCAL Hit Energy; Energy [MeV]; Hits / 100 MeV", 100, 0, 10000 );
  m_hitETot = new TH1I( "hitETot", "FCAL Hit Total Energy; Energy [MeV]; Events / 100 MeV", 100, 0, 10000 );
  m_hitT = new TH1I( "hitT", "FCAL Hit Time; t [ns]; Hits / 4 ns", 100, 0, 400 );
  m_hitT0 = new TH1I( "hitT0", "FCAL Hit Energy Weighted Average Time; t_{0} [ns]; Events / 4 ns", 100, 0, 400 );
  m_hitTmT0 = new TH1I( "hitTmT0", "FCAL Hit Local Time; t-t_{0} [ns]; Hits / 0.4 ns", 100, -20, 20 );
  m_hitE2D = new TH2F( "hitE2D", "FCAL Hit Average Energy [MeV]; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );
  m_hitTmT02D = new TH2F( "hitTmT02D", "FCAL Hit Local Time [ns]; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );
  m_hitTmT0Sq2D = new TH2F( "hitTmT0Sq2D", "FCAL Hit Local Time^{2} [ns^{2}]; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );
  m_hitOcc2D = new TH2F( "hitOcc2D", "FCAL Hit Occupancy; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5 );

  //
  // these with "clus" are after clusterization
  //
  
  m_clusN = new TH1I( "clusN", "FCAL Number of Clusters; Number of Clusters; Events", 40, -0.5, 39.5 );
  m_clusE = new TH1I( "clusE", "FCAL Cluster Energy; Energy [MeV]; Clusters / 50 MeV", 100, 0, 15000 );
  m_clusETot = new TH1I( "clusETot", "FCAL Cluster Total Energy; Energy [MeV]; Events / 100 MeV", 100, 0, 15000 );
  m_clusT = new TH1I( "clusT", "FCAL Cluster Time; t [ns]; Clusters / 4 ns", 100, 0, 400 );
  m_clusT0 = new TH1I( "clusT0", "FCAL Cluster Energy Weighted Average Time; t_{0} [ns]; Events / 4 ns]", 100, 0, 400 );
  m_clusTmT0 = new TH1I( "clusTmT0", "FCAL Cluster Local Time; t - t_{0} [ns]; Clusters / 0.8 ns", 100, -40, 40 );
  m_clusXYHigh = new TH2I( "clusXYHigh", "FCAL Cluster Positions (E > 200 MeV); x [cm]; y [cm]", 100, -150, 150, 100, -150, 150 );
  m_clusXYLow = new TH2I( "clusXYLow", "FCAL Cluster Positions (E < 200 MeV); x [cm]; y [cm]", 100, -150, 150, 100, -150, 150 );
  m_clusPhi = new TH1I( "clusPhi", "FCAL Cluster #phi; #phi [rad]; Clusters / 62.8 mrad", 100, -3.14, 3.14 );
  m_clusPhi->SetMinimum( 0 );
  m_clus2GMass = new TH1I( "clus2GMass", "FCAL 2 Cluster Invariant Mass E > 1 GeV; Invariant Mass [GeV]", 500, 0.0, 1.0 );

  //
  // these with "show" are after energy dependent non-linear correction
  //

  m_show2GMass = new TH1I( "show2GMass", "FCAL 2 Shower Invariant Mass E > 1 GeV; Invariant Mass [GeV]", 500, 0.0, 1.0 );
  m_showZvsE = new TH2I( "showZvsE", "FCAL z_{shower} vs. E_{shower}; E_{shower} [GeV]; z_{shower} [cm]", 120, 0.0, 6.0, 120, 620, 680 );
  m_showECorVsE = new TH2I( "showECorVsE", "FCAL E_{shower}/E_{cluster} vs. E_{shower}; E_{shower} [GeV]; E_{shower}/E_{cluster}", 200, 0.0, 6.0, 200, 0.0, 2.0 );
  m_showTsMTcVsZ = new TH2I( "showTsMTcVsZ", "FCAL t_{shower} - t_{cluster} vs. z_{shower}; z_{shower} cm; t_{shower} - t_{cluster} [ns]", 120, 620, 680, 120, -6, 6 );

  // back to main dir
  main->cd();

  // unlock
  japp->RootUnLock();

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_FCAL_online::brun(JEventLoop *eventLoop, int runnumber) {
  // This is called whenever the run number changes

  DGeometry *dgeom = NULL;
  DApplication *dapp = dynamic_cast< DApplication* >( eventLoop->GetJApplication() );
  if( dapp ) dgeom = dapp->GetDGeometry( runnumber );
   
  if( dgeom ){

    dgeom->GetTargetZ( m_targetZ );

  }
  else{

    cerr << "No geometry accessbile to FCAL_online monitoring plugin." << endl;
    return RESOURCE_UNAVAILABLE;
  }

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_FCAL_online::evnt(JEventLoop *eventLoop, int eventnumber) {

  vector< const DFCALGeometry* > geomVec;
  vector< const DFCALDigiHit*  > digiHits;
  vector< const DFCALHit*      > hits;
  vector< const DFCALCluster*  > clusterVec;
  vector< const DFCALShower*   > showerVec;
  eventLoop->Get( geomVec );
  eventLoop->Get( digiHits );
  eventLoop->Get( hits );
  if( hits.size() <= 500 ){  // only form clusters and showers if there aren't too many hits
    eventLoop->Get( clusterVec );
    eventLoop->Get( showerVec );
  }

  const DFCALGeometry& fcalGeom = *(geomVec[0]);

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //
  // extract and fill histograms for raw block level
  // hits -- this should be very fast
  //

  // loop through the hits and find the integral-weighted
  // average time for this event -- this is like a 
  // local t0 for the FCAL and will be useful
  // for looking at timing variances as there is no
  // absolute time in actual events
  
  double totIntegral = 0;
  double intWeightTime = 0;

  // while doing this, we find the Df250PulseIntegral and Df250PulsePedestal
  // objects (if any) associated with each hit and store them in a map.
  // This is done here since it is outside of the ROOT lock. It turns out
  // that the GetSingle calls use dynamic_cast and are therefore raher
  // expensive causing significant performance degredation if done inside
  // the lock. This will eventually be fixed in JANA, but for now we implement
  // this caching solution to make multi-threading efficient.  3/4/2015 DL
  map< const DFCALDigiHit*, pair<const Df250PulseIntegral*, const Df250PulsePedestal*> > pi_pp_cache;

  for( vector< const DFCALDigiHit* >::const_iterator dHitItr = digiHits.begin();
       dHitItr != digiHits.end(); ++dHitItr ){
  
    totIntegral += (**dHitItr).pulse_integral;
    intWeightTime += (**dHitItr).pulse_integral * (**dHitItr).pulse_time;

    // fetch lower level FADC data
    const Df250PulseIntegral* pulseInt = NULL;
    const Df250PulsePedestal* pulsePed = NULL;

    const DFCALDigiHit& dHit = (**dHitItr);
    dHit.GetSingle( pulseInt );
    if( pulseInt ) pulseInt->GetSingle( pulsePed );
	 
	 pi_pp_cache[&dHit] = pair<const Df250PulseIntegral*, const Df250PulsePedestal*>(pulseInt, pulsePed);
  }

  intWeightTime /= totIntegral;

  int nOverflow = 0;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // first find and energy weighted average time
  // for the hits in this event -- we'll use this
  // as a t0 later

  double hitETot = 0;
  double hitEwtT = 0;

  for( vector< const DFCALHit* >::const_iterator hit_itr = hits.begin();
       hit_itr != hits.end(); ++hit_itr ){

    hitETot += (**hit_itr).E;
    hitEwtT += (**hit_itr).E * (**hit_itr).t;
  }
  
  hitEwtT /= hitETot;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  double clusETot = 0;
  double clusEwtT = 0;

  for( vector< const DFCALCluster*>::const_iterator clusItr = clusterVec.begin();
       clusItr != clusterVec.end(); ++clusItr ){

    clusETot += (**clusItr).getEnergy();
    clusEwtT += (**clusItr).getEnergy() * (**clusItr).getTime();
  }

  clusEwtT /= clusETot;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  japp->RootWriteLock();

  if( digiHits.size() > 0 )
	  fcal_num_events->Fill(1);

  m_digT0->Fill( intWeightTime );
  m_digN->Fill( digiHits.size() );
  for( vector< const DFCALDigiHit* >::const_iterator dHitItr = digiHits.begin();
       dHitItr != digiHits.end(); ++dHitItr ){

    const DFCALDigiHit& dHit = (**dHitItr);

    int chan = fcalGeom.channel( dHit.row, dHit.column );

    m_digOcc2D->Fill( dHit.column, dHit.row );

    // fetch lower level FADC data
    pair<const Df250PulseIntegral*, const Df250PulsePedestal*> &pulseX = pi_pp_cache[&dHit];
    const Df250PulseIntegral* pulseInt = pulseX.first;
    const Df250PulsePedestal* pulsePed = pulseX.second;

    // dHit.GetSingle( pulseInt );
    // if( pulseInt ) pulseInt->GetSingle( pulsePed );

    if( pulseInt && pulsePed ){
      
      double peakV = pulsePed->pulse_peak * 2.0 / 4096;
      
      m_digPeakV->Fill( peakV );
      m_digPeakV2D->Fill( dHit.column, dHit.row, peakV );

      // overflows should be tagged at a lower level this
      // is an imperfect way of trying to detect them as
      // "pulse_peak" is defined as the first sample for
      // which the subsequent sample is a lower amplitude
      if( pulsePed->pulse_peak == 4096 ) ++nOverflow;

      m_digIntVsPeak->Fill( pulsePed->pulse_peak, pulseInt->integral );

      double integral = (double)( pulseInt->integral -  
				  ( pulseInt->pedestal * 
				    pulseInt->nsamples_integral ) );
      double peak = (double)( pulsePed->pulse_peak - pulseInt->pedestal );

      if( pulsePed->pulse_peak > 300 )
	m_digIntToPeak->Fill( integral / peak );
    }

    // pulse_time is a 15-bit int - the lower six bits are
    // a precision time determined by a separate algorithm
    // from the upper 9 bits:  track each separately

    int coarseT = dHit.pulse_time >> 6;
    int preciseT = dHit.pulse_time & 0x3F;

    m_digInt->Fill( dHit.pulse_integral );
    m_digT->Fill( dHit.pulse_time );
    if( digiHits.size() > 0 )
      m_digTmT0->Fill( dHit.pulse_time - intWeightTime );
    m_digTmT02D->Fill( dHit.column, dHit.row, dHit.pulse_time - intWeightTime );
    m_digCoarseT->Fill( coarseT );
    m_digCoarseTChan->Fill( chan, coarseT );
    m_digPreciseT->Fill( preciseT );
    m_digPreciseTChan->Fill( chan, preciseT );
    m_digPed->Fill( dHit.pedestal );
    m_digPedChan->Fill( chan, dHit.pedestal );
    m_digPed2D->Fill( dHit.column, dHit.row, dHit.pedestal );
    m_digPedSq2D->Fill( dHit.column, dHit.row, dHit.pedestal*dHit.pedestal );
    m_digQual->Fill( dHit.QF );
  }
  m_digNOver->Fill( nOverflow );

//  japp->RootUnLock();

  // end of raw hit filling

  // 
  // extract and fill histograms for calibrated block
  // level hits -- this should also be fast
  // 

  
  // now fill histograms

//  japp->RootWriteLock();

  m_hitETot->Fill( hitETot*k_to_MeV );
  m_hitT0->Fill( hitEwtT );

  m_hitN->Fill( hits.size() );
  for( vector< const DFCALHit* >::const_iterator hit_itr = hits.begin();
       hit_itr != hits.end(); ++hit_itr ){

    const DFCALHit& hit = (**hit_itr);

    double locTime = ( hit.t - hitEwtT )*k_to_nsec;

    m_hitE->Fill( hit.E*k_to_MeV );
    m_hitT->Fill( hit.t*k_to_nsec );

    if( hits.size() > 1 )
      m_hitTmT0->Fill( locTime );

    m_hitE2D->Fill( hit.column, hit.row, hit.E*k_to_MeV );
    m_hitTmT02D->Fill( hit.column, hit.row, locTime );
    m_hitTmT0Sq2D->Fill( hit.column, hit.row, locTime*locTime );
    m_hitOcc2D->Fill( hit.column, hit.row );
  }

//  japp->RootUnLock();

  // end calibrated block-level filling

  // for events with a lot of hits -- stop now
  if( hits.size() > 500 ){
    japp->RootUnLock();
    return NOERROR;
  }

  //
  // if there are a small number of hits go ahead
  // and run the clusterizer and make a few plots
  // utilizing the list of clusters and showers
  //




//  japp->RootWriteLock();
  
  m_clusN->Fill( clusterVec.size() );

  if( clusterVec.size() > 0 ){

    m_clusT0->Fill( clusEwtT );
    m_clusETot->Fill( clusETot * k_to_MeV );
  }

  for( vector< const DFCALCluster*>::const_iterator clusItr = clusterVec.begin();
       clusItr != clusterVec.end(); ++clusItr ){

    const DFCALCluster& cluster = (**clusItr);

    double clusX = cluster.getCentroid().X();
    double clusY = cluster.getCentroid().Y();

    double clusR = sqrt( clusX * clusX + clusY * clusY );

    DVector3 clus1Mom = cluster.getCentroid();
    clus1Mom.SetMag( cluster.getEnergy() );
    
    m_clusPhi->Fill( clus1Mom.Phi() );
    m_clusE->Fill( cluster.getEnergy() * k_to_MeV );
    m_clusT->Fill( cluster.getTime() );

    if( clusterVec.size() > 1 )
      m_clusTmT0->Fill( cluster.getTime() - clusEwtT );

    if( cluster.getEnergy() > 200*k_MeV ){

      m_clusXYHigh->Fill( cluster.getCentroid().X(),
			  cluster.getCentroid().Y() );
    }
    else{

      m_clusXYLow->Fill( cluster.getCentroid().X(),
			 cluster.getCentroid().Y() );
    }

    if( ( cluster.getEnergy() < 1*k_GeV ) ||
	( clusR < 20*k_cm ) ||
	( cluster.GetNHits() < 2 )
	) continue;

    for( vector< const DFCALCluster*>::const_iterator clus2Itr = clusItr + 1;
	 clus2Itr != clusterVec.end(); ++clus2Itr ){

      const DFCALCluster& cluster2 = (**clus2Itr);

      double clus2X = cluster2.getCentroid().X();
      double clus2Y = cluster2.getCentroid().Y();

      double clus2R = sqrt( clus2X * clus2X + clus2Y * clus2Y );
      
      // only make the 2-cluster invariant mass plot for
      // cases where both clusters are greater than 200 MeV
      // in energy

      if( ( cluster2.getEnergy() < 1*k_GeV ) ||
	  ( clus2R < 20*k_cm ) ||
	  ( cluster2.GetNHits() < 2 )
	  ) continue;

      if( fabs( cluster.getTime() - cluster2.getTime() ) > 10*k_nsec ) continue;
      
      DVector3 clus2Mom = cluster2.getCentroid();
      clus2Mom.SetMag( cluster2.getEnergy() );

      DLorentzVector gam1( clus1Mom, clus1Mom.Mag() );
      DLorentzVector gam2( clus2Mom, clus2Mom.Mag() );

      m_clus2GMass->Fill( ( gam1 + gam2 ).M() );
    }
  }
//  japp->RootUnLock();

  // end cluster lever filling

  // now do shower filling -- these are clusters that have
  // had timing, position, and nonlinear energy corrections
  // applied on a cluster level


//  japp->RootWriteLock();

  for( vector< const DFCALShower* >::const_iterator showItr = showerVec.begin();
       showItr != showerVec.end(); ++showItr ){

    const DFCALShower& show = (**showItr);

    const DFCALCluster* clus = NULL;
    show.GetSingle( clus );
    if( clus == NULL ) continue;

    DVector3 clusMom = clus->getCentroid();
    clusMom.SetZ( clusMom.Z() - m_targetZ );
    clusMom.SetMag( clus->getEnergy() );
 
    DVector3 showMom = show.getPosition();
    showMom.SetZ( showMom.Z() - m_targetZ );
    showMom.SetMag( show.getEnergy() );
    
    m_showZvsE->Fill( show.getEnergy(), show.getPosition().Z() );
    m_showECorVsE->Fill( show.getEnergy(), show.getEnergy() / clus->getEnergy() );
    m_showTsMTcVsZ->Fill( show.getPosition().Z(), show.getTime() - clus->getTime() );

    if( ( show.getEnergy() < 1*k_GeV ) ||
	( show.getPosition().Pt() < 20*k_cm ) ) continue;

    for( vector< const DFCALShower* >::const_iterator show2Itr = showItr + 1;
	 show2Itr != showerVec.end(); ++show2Itr ){

      const DFCALShower& show2 = (**show2Itr);
      
      // only make the 2-shower invariant mass plot for
      // cases where both showers are greater than 200 MeV
      // in energy

      if( ( show2.getEnergy() < 1*k_GeV ) ||
	  ( show2.getPosition().Pt() < 20*k_cm ) ) continue; 

      DVector3 show2Mom = show2.getPosition();
      show2Mom.SetZ( show2Mom.Z() - m_targetZ );
      show2Mom.SetMag( show2.getEnergy() );

      DLorentzVector gam1( showMom, showMom.Mag() );
      DLorentzVector gam2( show2Mom, show2Mom.Mag() );

      m_show2GMass->Fill( ( gam1 + gam2 ).M() );
    }
    

  }
  japp->RootUnLock();
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_FCAL_online::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_FCAL_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
