// $Id$
//
//    File: DFCALHit_factory.cc
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <cassert>
using namespace std;

#include "FCAL/DFCALHit_factory.h"
#include "FCAL/DFCALHit.h"
#include "FCAL/DFCALMCResponse.h"
#include "FCAL/DFCALGeometry.h"

//------------------
// evnt
//------------------
jerror_t DFCALHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	
	assert( _data.size() == 0 );
	
	// extract the FCAL Geometry
	vector<const DFCALGeometry*> fcalGeomVect;
	loop->Get( fcalGeomVect );
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);

    vector<const DFCALMCResponse*> responseVect;
    loop->Get( responseVect );
    
    for( vector<const DFCALMCResponse*>::const_iterator resp = responseVect.begin();
         resp != responseVect.end();
         ++resp ){
     
        DFCALHit* hit = new DFCALHit();
        
        hit->id = (**resp).id;
        hit->E = (**resp).E();
        hit->t = (**resp).t();
        
        TVector2 pos = fcalGeom.positionOnFace( (**resp).channel() );
        hit->x = pos.X();
        hit->y = pos.Y();
        
        _data.push_back( hit );
    }
    
	return NOERROR;
}
    
