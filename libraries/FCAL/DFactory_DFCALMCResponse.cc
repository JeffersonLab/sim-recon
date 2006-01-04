// $Id$
//
//    File: DFactory_DFCALMCResponse.cc
// Created: Wed Aug 24 10:05:31 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#include <cassert>

#include "DFactory_DFCALMCResponse.h"
#include "DFCALTruthShower.h"
#include "DFCALGeometry.h"

//------------------
// evnt
//------------------
derror_t DFactory_DFCALMCResponse::evnt(DEventLoop *loop, int eventnumber)
{
	
	assert( _data.size() == 0 );
	
	// extract the FCAL Geometry
	vector<const DFCALGeometry*> fcalGeomVect;
	loop->Get( fcalGeomVect );
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);
	
	// extract the HDDM objects
	vector<const DFCALTruthShower*> truthShowerVect;
	loop->Get( truthShowerVect );
	
	for( vector<const DFCALTruthShower*>::const_iterator 
		 truthSh = truthShowerVect.begin();
		 truthSh != truthShowerVect.end();
		 ++truthSh ){
		
		// loop over each HDDM shower and make the FCALMCResponse
		// objects for showers that hit real FCAL blocks
		
		int row = fcalGeom.row( (**truthSh).y() );
		int col = fcalGeom.column( (**truthSh).x() );
		
		if( fcalGeom.isBlockActive( row, col ) ){
		
			// can do something fancier later -- like smear E and t
			
			_data.push_back
			  ( new DFCALMCResponse( (**truthSh).id,
									 fcalGeom.channel( row, col ),
									 (**truthSh).E(),
									 (**truthSh).t() ) );
		}
	}
	
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DFCALMCResponse::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	//		printheader("row:    x:     y:");
	//
	// 	for(int i=0; i<_data.size(); i++){
	//			DFCALMCResponse *myDFCALMCResponse = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDFCALMCResponse->x);
	//			printcol("%3.2f",	myDFCALMCResponse->y);
	//			printrow();
	//		}
	//
	return _table;

}

