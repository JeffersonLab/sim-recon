// $Id$
//
//    File: DFCALHit_factory.cc
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <cassert>

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
    
//------------------
// toString
//------------------
const string DFCALHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("hit:  column:   row:   x(cm):   y(cm):   E(MeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALHit *fcalhit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%3.1f", fcalhit->x);
		printcol("%3.1f", fcalhit->y);
		printcol("%2.3f", fcalhit->E*1000.0);
		printcol("%4.0f", fcalhit->t);
		printrow();
	}
	
	return _table;
}
