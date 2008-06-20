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
    
    // set cell threshold at 2.5 MeV
    m_cellThreshold = 2.5 * k_MeV;

    // set the sampling smearing coefficients:
    // (from GlueX-doc 827 v3 Figure 13 )
    m_samplingCoefA = 0.042;
    m_samplingCoefB = 0.013;
    
    gPARMS->SetDefaultParameter( "BCALRESPONSE:CELL_THRESHOLD",  m_cellThreshold );
    gPARMS->SetDefaultParameter( "BCALRESPONSE:SAMPLING_COEF_A", m_samplingCoefA );
    gPARMS->SetDefaultParameter( "BCALRESPONSE:SAMPLING_COEF_B", m_samplingCoefB );
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

        float smearedE = samplingSmear( hddmhit->E );
        
        float upDist = ( bcalGeom.BCALFIBERLENGTH / 2 ) + hddmhit->zLocal;
        float downDist = ( bcalGeom.BCALFIBERLENGTH / 2 ) - hddmhit->zLocal;
        
        // sampling fluctuations are correlated between ends
        float upEnergy = smearedE * exp( -upDist / bcalGeom.ATTEN_LENGTH );
        float downEnergy = smearedE * exp( -downDist / bcalGeom.ATTEN_LENGTH );

        float upTime = hddmhit->t + upDist / bcalGeom.C_EFFECTIVE;
        float downTime = hddmhit->t + downDist / bcalGeom.C_EFFECTIVE;
        
        if( upEnergy > m_cellThreshold ){
        
            DBCALMCResponse *response = new DBCALMCResponse;
        
            response->module =hddmhit->module;
            response->layer = hddmhit->layer;
            response->sector = hddmhit->sector;
            response->E = upEnergy;
            response->t = upTime;
            response->end = DBCALGeometry::kUpstream;

            response->cellId = 
                DBCALGeometry::cellId( hddmhit->module, hddmhit->layer, hddmhit->sector );

            _data.push_back(response);
        }
        
        if( downEnergy > m_cellThreshold ){
            
            DBCALMCResponse *response = new DBCALMCResponse;
            
            response->module =hddmhit->module;
            response->layer = hddmhit->layer;
            response->sector = hddmhit->sector;
            response->E = downEnergy;
            response->t = downTime;
            response->end = DBCALGeometry::kDownstream;
            
            response->cellId = 
                DBCALGeometry::cellId( hddmhit->module, hddmhit->layer, hddmhit->sector );
            
            _data.push_back(response);
        }
    }
    
	return NOERROR;
}

float
DBCALMCResponse_factory::samplingSmear( float E )
{
    double sigma = m_samplingCoefA / sqrt( E ) + m_samplingCoefB;
    
    return( E * m_randomGen.Gaus( 1., sigma ) );
}

