// $Id$
//
//    File: DBCALMCResponse_factory.cc
// Created: Thu Nov 17 09:56:05 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#include <cassert>

#include "BCAL/DBCALMCResponse_factory.h"
#include "BCAL/DBCALGeometry.h"
#include "BCAL/DHDDMBCALHit.h"

#include "DRandom.h"
#include "units.h"

DBCALMCResponse_factory::DBCALMCResponse_factory() : 
m_randomGen()
{
    
    // setup response parameters
    
    // set cell threshold at 10 MeV
    //m_cellThreshold = 2.5 * k_MeV;
    m_cellThreshold = 10 * k_MeV;


    // set the sampling smearing coefficients:
    // (from http://www.jlab.org/Hall-D/software/wiki/index.php/Image:ReSAM_samfrac_1.gif )
    m_samplingCoefA = 0.041;
    m_samplingCoefB = 0.09;
}

//------------------
// evnt
//------------------
jerror_t DBCALMCResponse_factory::evnt(JEventLoop *loop, int eventnumber)
{

    vector<const DBCALGeometry*> bcalGeomVec;
    eventLoop->Get(bcalGeomVec);
    const DBCALGeometry& bcalGeom = *(bcalGeomVec.at( 0 ));
    
    vector<const DHDDMBCALHit*> hddmhits;
    eventLoop->Get(hddmhits);
    
    for (unsigned int i = 0; i < hddmhits.size(); i++) {
        
        const DHDDMBCALHit *hddmhit = hddmhits[i];
        DBCALMCResponse *response = new DBCALMCResponse;
        
        response->module =hddmhit->module;
        response->layer = hddmhit->layer;
        response->sector = hddmhit->sector;
        response->E = samplingSmear( hddmhit->E );
        response->t = hddmhit->t;

        response->cellId = 
            DBCALGeometry::cellId( hddmhit->module, hddmhit->layer, hddmhit->sector );

        float upDist =  ( bcalGeom.BCALFIBERLENGTH / 2 ) + hddmhit->zLocal;
        float downDist = ( bcalGeom.BCALFIBERLENGTH / 2 ) - hddmhit->zLocal;
        
        
        // make sure we exceed the effective electronics threshold
        if( response->E > m_cellThreshold ){
            
            _data.push_back(response);
        }
    }
    
	return NOERROR;
}

float
DBCALMCResponse_factory::samplingSmear( float E )
{
 
    double sigma = sqrt( ( m_samplingCoefA * m_samplingCoefA / E ) +
                         ( m_samplingCoefB * m_samplingCoefB ) );
    
    return( E * m_randomGen.Gaus( 1., sigma ) );
}

//------------------
// toString
//------------------
const string DBCALMCResponse_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The JFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	//		printheader("row:    x:     y:");
	//
	// 	for(int i=0; i<_data.size(); i++){
	//			DBCALMCResponse *myDBCALMCResponse = _data[i];
	//
	//			printnewrow();
	//			printcol("%d",	i);
	//			printcol("%1.3f",	myDBCALMCResponse->x);
	//			printcol("%3.2f",	myDBCALMCResponse->y);
	//			printrow();
	//		}
	//
	printheader("row:   module:  layer:  sector:         end:     E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DBCALMCResponse *s = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%d",	s->module);
		printcol("%d",	s->layer);
		printcol("%d",	s->sector);
		printcol(s->end==0? "upstream":"downstream");
		printcol("%2.3f",	s->E);
		printcol("%2.3f",	s->t);
		printrow();
	}


	return _table;

}
