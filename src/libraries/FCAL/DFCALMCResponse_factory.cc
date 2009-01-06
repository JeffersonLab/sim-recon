// $Id$
//
//    File: DFCALMCResponse_factory.cc
// Created: Wed Aug 24 10:05:31 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#include <cassert>
using namespace std;

#include "FCAL/DFCALMCResponse_factory.h"
#include "FCAL/DMCFCALHit.h"
#include "FCAL/DFCALGeometry.h"

DFCALMCResponse_factory::DFCALMCResponse_factory() :
m_randomGen()
{

    // setup response parameters

    // set block  threshold at 18 MeV 
    // this is ~ equivalent to 30 MeV incident energy for 166 cm attenuation length
    m_blockThreshold = 20.0 * k_MeV;

    // set the photon-statistics factor for smearing hit energy, from Criss's MC
    m_photStatCoef = 0.035;

    gPARMS->SetDefaultParameter( "FCALRESPONSE:BLOCK_THRESHOLD",  m_blockThreshold );
    gPARMS->SetDefaultParameter( "FCALRESPONSE:PHOT_STAT_C", m_photStatCoef );
}


//------------------
// evnt
//------------------
jerror_t DFCALMCResponse_factory::evnt(JEventLoop *loop, int eventnumber)
{
	
	assert( _data.size() == 0 );
	
	// extract the FCAL Geometry
	vector<const DFCALGeometry*> fcalGeomVect;
	loop->Get( fcalGeomVect );
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);
	
	// extract the HDDM objects
	vector<const DMCFCALHit*> hitVect;
	loop->Get( hitVect );
	
	for( vector<const DMCFCALHit*>::const_iterator 
		 hit = hitVect.begin();
		 hit != hitVect.end();
		 ++hit ){
		
		// loop over each HDDM shower and make the FCALMCResponse
		// objects for showers that hit real FCAL blocks
		
		int row = (**hit).row;
		int col = (**hit).column;
		
		if( fcalGeom.isBlockActive( row, col ) ){
		
			// can do something fancier later -- like smear E and t
			// or implement Cherenkov photon effects

                        double sigma = m_photStatCoef / sqrt( (**hit).E ) ;
                	float smearE = (**hit).E * m_randomGen.Gaus( 1., sigma );;
                       
                        if ( smearE > m_blockThreshold ) {  
				_data.push_back
			  	( new DFCALMCResponse( (**hit).id,
									 fcalGeom.channel( row, col ),
									 smearE, // (**hit).E,
										 (**hit).t ) );
			}
		}
	}
	
	return NOERROR;
}

