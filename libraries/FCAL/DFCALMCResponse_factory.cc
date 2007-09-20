// $Id$
//
//    File: DFCALMCResponse_factory.cc
// Created: Wed Aug 24 10:05:31 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#include <cassert>

#include "FCAL/DFCALMCResponse_factory.h"
#include "FCAL/DMCFCALHit.h"
#include "FCAL/DFCALGeometry.h"

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
            
			_data.push_back
			  ( new DFCALMCResponse( (**hit).id,
									 fcalGeom.channel( row, col ),
									 (**hit).E,
									 (**hit).t ) );
		}
	}
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFCALMCResponse_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("id:      chanel: E(MeV): t(ns):");

	for(unsigned int i=0; i<_data.size(); i++){
		DFCALMCResponse *mchit = _data[i];

		printnewrow();
		printcol("%d",	mchit->id );
		printcol("%3.1f",	mchit->channel() );
		printcol("%3.1f",	mchit->E() );
		printcol("%3.1f",	mchit->t() );
		printrow();
	}

	return _table;

}

