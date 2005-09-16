// $Id$
//
//    File: DFactory_DFCALHit.cc
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <cassert>

#include "DFactory_DFCALHit.h"
#include "DFCALHit.h"
#include "DFCALMCResponse.h"
#include "DFCALGeometry.h"

//------------------
// evnt
//------------------
derror_t DFactory_DFCALHit::evnt(DEventLoop *eventLoop, int eventnumber)
{
	assert( _data.size() == 0 );
	
	vector<const DFCALMCResponse*> responseVect;
	eventLoop->Get( responseVect );
		
	// extract the FCAL Geometry
	vector<const DFCALGeometry*> fcalGeomVect;
	eventLoop->Get( fcalGeomVect );
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);
		
	for( vector<const DFCALMCResponse*>::const_iterator 
		 aResponse = responseVect.begin();
		 aResponse != responseVect.end();
		 ++aResponse ){
		
		// temporary for now:
	
		float x = fcalGeom.positionOnFace( (**aResponse).channel() ).X();
		float y = fcalGeom.positionOnFace( (**aResponse).channel() ).Y();
		
		_data.push_back( new DFCALHit( (**aResponse).id,
									   x, y,
									   (**aResponse).E(),
									   (**aResponse).t() ) );
	}
	
	return NOERROR;
}


//------------------
// toString
//------------------
const string DFactory_DFCALHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   x(cm):   y(cm):   E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALHit *fcalhit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%3.1f", fcalhit->x);
		printcol("%3.1f", fcalhit->y);
		printcol("%2.3f", fcalhit->E);
		printcol("%4.0f", fcalhit->t);
		printrow();
	}
	
	return _table;
}
