// $Id$
//
//    File: DBCALMCResponse_factory.cc
// Created: Thu Nov 17 09:56:05 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#include <cassert>
#include <iostream>

#include "BCAL/DBCALMCResponse_factory.h"
#include "BCAL/DBCALGeometry.h"
#include "BCAL/DHDDMBCALHit.h"

#include "DRandom.h"
#include "units.h"

using namespace std;

DBCALMCResponse_factory::DBCALMCResponse_factory() : 
m_randomGen()
{
    
  // setup response parameters
  m_darkRate_GHz         = 0.041;
  m_xTalk_fract          = 0.03;
  m_intWindow_ns         = 100;
  m_devicePDE            = 0.12;
  m_sampling_fract       = 0.15;
  m_maxOccupancy_fract   = 0.05;

  // GX-doc 1069, Table 1 -- try to extract back to
  // photons per side per MeV in fiber
  // 4.6 / PDE / attentuation  (meaurements performed in center)
  // 75 = 4.6  / 0.12 / exp( -200 / 300 )
  m_photonsPerSidePerMeVInFiber = 75;

  // set the sampling smearing coefficients:
  // (from GlueX-doc 827 v3 Figure 13 )
  m_samplingCoefA = 0.042;
  m_samplingCoefB = 0.013;

  // time smearing comes from beam test time difference resolution with
  // beam incident on the center of the module
  m_timediffCoefA = 0.07 * sqrt( 2 );
  // no floor term, but leave the option here:
  m_timediffCoefB = 0.0 * sqrt( 2 );
    
  // set this low for now -- needs more thought later
  m_cellThresholdOuter = 1 * k_MeV;
  
  gPARMS->SetDefaultParameter( "BCALRESPONSE:CELL_THRESHOLD_OUTER",  m_cellThresholdOuter );

  gPARMS->SetDefaultParameter( "BCALRESPONSE:DARK_RATE_GHZ", m_darkRate_GHz );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:CROSS_TALK_PROB", m_xTalk_fract );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:FADC_WINDOW_NS", m_intWindow_ns );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:DEVICE_PDE", m_devicePDE );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:SAMPLING_FRACTION", m_sampling_fract );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:OCCUPANCY_FRACTION_LIMIT", 
			       m_maxOccupancy_fract );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:PHOTONS_PER_SIDE_PER_MEV_IN_FIBER", 
			       m_photonsPerSidePerMeVInFiber );

  gPARMS->SetDefaultParameter( "BCALRESPONSE:SAMPLING_COEF_A", m_samplingCoefA );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:SAMPLING_COEF_B", m_samplingCoefB );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:TIMESMEAR_COEF_A", m_timediffCoefA );
  gPARMS->SetDefaultParameter( "BCALRESPONSE:TIMESMEAR_COEF_B", m_timediffCoefB );

}

jerror_t
DBCALMCResponse_factory::init( void )
{

  // this the average number of (assumed single PE) 
  // dark pulses in a window

  double nAvg = m_darkRate_GHz * m_intWindow_ns;

  // now we need to find n such that 
  // sum_i=1^n P(n) > ( 1 - maxOccupancy )
  // P(n) is the probability to have n pulses
  //
  // P(n) is really a convolution since we have:
  // P(n) = P( n_d + n_x ) = 
  //    P( n_d; nAvg ) * P( n_x; n_d * x )
  // 
  // n_d number of dark pulses
  // x is the cross talk rate

  // numerically build int_0^n P(n)

  // in practice the cutoff is going to be < 100 PE
  // we can build an accurate pdf up to that
  double pdf[100];
  for( int i = 0; i < 100; ++i ){ pdf[i] = 0; }

  // probability for zero is easy:
  double darkTerm = exp( -nAvg );
  pdf[0] = darkTerm;

  for( int n_d = 1; n_d < 100; ++n_d ){

    darkTerm *= ( nAvg / n_d );

    double xTalkAvg = n_d * m_xTalk_fract;

    // probability for zero x-talk pulses
    double xTerm = exp( -xTalkAvg );
    pdf[n_d] += ( xTerm * darkTerm );

    // now include probability for additional
    // cross talk pulses
    for( int n_x = 1; n_x + n_d < 100; ++n_x ){

      xTerm *= ( xTalkAvg / n_x );

      pdf[n_d+n_x] += ( xTerm * darkTerm );
    }
  }

  double integral = 0;
  int nPEThresh = 0;
  while( integral < ( 1 - m_maxOccupancy_fract ) ){

    // post increment includes zero and requires
    // one more PE than what breaks the loop
    integral += pdf[nPEThresh];
    ++nPEThresh;
  }

  // now get the photon theshold
  double photonThresh = nPEThresh / m_devicePDE;

  // now convert this into a energy threshold
  m_cellThresholdInner = 
    ( photonThresh / m_photonsPerSidePerMeVInFiber ) / 
    m_sampling_fract * k_MeV; 

  cout << "BCAL inner cell threshold:  " 
       << m_cellThresholdInner * 1000 << " MeV, "  
       << nPEThresh << " P.E." << endl;

  return NOERROR;
}


//------------------
// evnt
//------------------
jerror_t DBCALMCResponse_factory::evnt(JEventLoop *loop, int eventnumber)
{
 
    vector<const DBCALGeometry*> bcalGeomVec;
    loop->Get(bcalGeomVec);
    const DBCALGeometry& bcalGeom = *(bcalGeomVec.at( 0 ));
    
    map< int, pair< int, int > > darkHits;
    map< int, pair< int, int > >::iterator darkHitItr;

    for( int m = 1; m <= bcalGeom.NBCALMODS; ++m ){
      for( int l = 1; l <= bcalGeom.NBCALLAYS1; ++l ){
	for( int s = 1; s <= bcalGeom.NBCALSECS1; ++s ){

	  pair< int, int > nHits( getDarkHits(), getDarkHits() );

	  darkHits[DBCALGeometry::cellId( m, l, s )] = nHits;
	}
      }
    }

    double mevPerPE = 1 / 
      ( m_photonsPerSidePerMeVInFiber * m_devicePDE *
	m_sampling_fract );

    vector<const DHDDMBCALHit*> hddmhits;
    loop->Get(hddmhits);
    
    for (unsigned int i = 0; i < hddmhits.size(); i++) {

        const DHDDMBCALHit *hddmhit = hddmhits[i];

        float smearedE = samplingSmear( hddmhit->E );

        float upDist = ( bcalGeom.BCALFIBERLENGTH / 2 ) + hddmhit->zLocal;
        float downDist = ( bcalGeom.BCALFIBERLENGTH / 2 ) - hddmhit->zLocal;
        
        // sampling fluctuations are correlated between ends
        float upEnergy = smearedE * exp( -upDist / bcalGeom.ATTEN_LENGTH );
        float downEnergy = smearedE * exp( -downDist / bcalGeom.ATTEN_LENGTH );

	// independently smear time for both ends -- time smearing 
	// parameters come from data taken with beam at the center of 
	// the module so there is an implcit exp( ( -L / 2 ) / lambda ) 
	// that needs to be canceled out since we are working
	// at this stage with attenuated energies 
        float smearedtUp = 
	  timeSmear( hddmhit->t, 
		     upEnergy * exp( ( bcalGeom.BCALFIBERLENGTH / 2 ) / 
				     bcalGeom.ATTEN_LENGTH ) );
	
        float smearedtDown = 
	  timeSmear( hddmhit->t, 
		     downEnergy * exp( ( bcalGeom.BCALFIBERLENGTH / 2 ) / 
				       bcalGeom.ATTEN_LENGTH ) );

	// now offset times for propagation distance
	float upTime = smearedtUp + upDist / bcalGeom.C_EFFECTIVE;
	float downTime = smearedtDown + downDist / bcalGeom.C_EFFECTIVE;
	
	// hdgeat numbers layers starting with 1
	float cellThreshold = ( hddmhit->layer <= bcalGeom.NBCALLAYS1 ?
				m_cellThresholdInner : m_cellThresholdOuter );

	int cell = DBCALGeometry::cellId( hddmhit->module, hddmhit->layer, 
					  hddmhit->sector );

	darkHitItr = darkHits.find( cell );
	if( darkHitItr != darkHits.end() ){

	  upEnergy += ( darkHitItr->second.first * mevPerPE * k_MeV );
	  downEnergy += ( darkHitItr->second.second * mevPerPE * k_MeV );

	  // now delete this from the map so we don't create
	  // additional hits later

	  darkHits.erase( darkHitItr );
	}
	
        if( upEnergy > cellThreshold ){
        
            DBCALMCResponse *response = new DBCALMCResponse;
        
            response->module =hddmhit->module;
            response->layer = hddmhit->layer;
            response->sector = hddmhit->sector;
            response->E = upEnergy;
            response->t = upTime;
            response->end = DBCALGeometry::kUpstream;

            response->cellId = cell;

            _data.push_back(response);
        }
        
        if( downEnergy > cellThreshold ){
            
            DBCALMCResponse *response = new DBCALMCResponse;
            
            response->module =hddmhit->module;
            response->layer = hddmhit->layer;
            response->sector = hddmhit->sector;
            response->E = downEnergy;
            response->t = downTime;
            response->end = DBCALGeometry::kDownstream;
	    response->cellId = cell;
            
            _data.push_back(response);
        }
    }

    
    // make hits for remainig dark pulses...
    for( darkHitItr = darkHits.begin();
	 darkHitItr != darkHits.end();
	 ++darkHitItr ){
      
      int cell = darkHitItr->first;
      double upEnergy = darkHitItr->second.first * mevPerPE * k_MeV;
      double downEnergy = darkHitItr->second.second * mevPerPE * k_MeV;
      
      if( upEnergy > m_cellThresholdInner ){
        
	DBCALMCResponse *response = new DBCALMCResponse;
        
	response->module = DBCALGeometry::module( cell );
	response->layer = DBCALGeometry::layer( cell );
	response->sector = DBCALGeometry::sector( cell );
	response->E = upEnergy;
	response->t = m_randomGen.Uniform( -0.25 * m_intWindow_ns,
					    0.75 * m_intWindow_ns ) * k_nsec;
	response->end = DBCALGeometry::kUpstream;
	    
	response->cellId = cell;

	_data.push_back(response);
      }
        
      if( downEnergy > m_cellThresholdInner ){
            
	DBCALMCResponse *response = new DBCALMCResponse;
        
	response->module = DBCALGeometry::module( cell );
	response->layer = DBCALGeometry::layer( cell );
	response->sector = DBCALGeometry::sector( cell );
	response->E = downEnergy;
	response->t = m_randomGen.Uniform( -0.25 * m_intWindow_ns,
					    0.75 * m_intWindow_ns ) * k_nsec;
	response->end = DBCALGeometry::kDownstream;
	response->cellId = cell;
        
	_data.push_back(response);
      }
    }
    
    return NOERROR;
}

float
DBCALMCResponse_factory::samplingSmear( float E )
{
    double sigmaSamp = m_samplingCoefA / sqrt( E ) + m_samplingCoefB;
    
    return( E * m_randomGen.Gaus( 1., sigmaSamp ) );
}

float
DBCALMCResponse_factory::timeSmear( float t, float E )
{

  double sigmaT = m_timediffCoefA / sqrt( E ) + m_timediffCoefB;

  return( t + m_randomGen.Gaus( 0., sigmaT ) );
}

int
DBCALMCResponse_factory::getDarkHits()
{

  int darkPulse = 
    m_randomGen.Poisson( m_darkRate_GHz * m_intWindow_ns );

  int xTalk = m_randomGen.Poisson( darkPulse * m_xTalk_fract );

  return( xTalk + darkPulse );
}
