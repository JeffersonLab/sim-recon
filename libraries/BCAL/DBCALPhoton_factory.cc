/*
 *  DBCALPhoton_factory.cc
 *  Hall D
 *
 *  Created by Matthew Shepherd on 7/23/07.
 *
 */

#include "BCAL/DBCALPhoton_factory.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALGeometry.h"

DBCALPhoton_factory::DBCALPhoton_factory()
{
    
    // energy calibration parameters
    m_scaleZ_p0 =  8.015E-01;
    m_scaleZ_p1 =  1.740E-03;
    m_scaleZ_p2 = -6.804E-06;
    m_scaleZ_p3 =  7.162E-09;
    
    m_nonlinZ_p0 =  4.441E-02;
    m_nonlinZ_p1 = -1.183E-04;
    m_nonlinZ_p2 =  3.482E-07;
    
}

//------------------
// brun
//------------------
jerror_t DBCALPhoton_factory::brun(JEventLoop *loop, int runnumber)
{
    
    vector<const DBCALGeometry*> bcalGeomVect;
    eventLoop->Get( bcalGeomVect );
    const DBCALGeometry& bcalGeom = *(bcalGeomVect[0]);

    m_bcalIR = bcalGeom.BCALINNERRAD;
    m_zTarget = 65;                    // global target position -- should come from database!
    
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALPhoton_factory::evnt(JEventLoop *loop, int eventnumber)
{
    
    vector< const DBCALShower* > showerVect;
    loop->Get( showerVect, "IU" );
    
    for( vector< const DBCALShower* >::iterator showItr = showerVect.begin();
         showItr != showerVect.end();
         ++showItr ){
                
        double xSh = (**showItr).x;
        double ySh = (**showItr).y;
        double zSh = (**showItr).z;
        
        // get z where shower enters BCAL (this corresponds to generated z
        // in tuning MC)
        double zEntry = zSh - ( ( zSh - m_zTarget ) * 
                                ( 1 - m_bcalIR / ( sqrt( xSh * xSh + ySh * ySh ) ) ) ); 
        
        // calibrate energy:
        // Energy calibration has a z dependence -- the
        // calibration comes from fitting E_rec / E_gen to scale * E_gen^nonlin
        // for slices of z.  These fit parameters (scale and nonlin) are then plotted 
        // as a function of z and fit.
        

        double nonlin = m_nonlinZ_p0 +
            m_nonlinZ_p1 * zEntry +
            m_nonlinZ_p2 * zEntry * zEntry;

        // for the scale, we extend the scale factor at the target as a constant behind
        // the target
        if( zEntry < m_zTarget ) zEntry = m_zTarget;
        
        double scale = m_scaleZ_p0 +
            m_scaleZ_p1 * zEntry +
            m_scaleZ_p2 * zEntry * zEntry +
            m_scaleZ_p3 * zEntry * zEntry * zEntry;

        
        // now turn E_rec into E_gen -->> E_gen = ( E_rec / scale ) ^ ( 1 / ( 1 + nonlin ) )
        double energy = pow( (**showItr).E / scale, 1 / ( 1 + nonlin ) );
        
        double pScale = energy / 
            sqrt( xSh * xSh + ySh * ySh + ( zSh - m_zTarget ) * ( zSh - m_zTarget ) );
        
        DBCALPhoton* photon = new DBCALPhoton( (**showItr).id );
        
        photon->setLorentzMomentum( DLorentzVector( xSh * pScale, 
                                                    ySh * pScale, 
                                                    ( zSh - m_zTarget ) * pScale, 
                                                    energy ) );
        
        photon->setShowerPosition( DVector3( xSh, ySh, zSh ) );
        
        photon->setFitLayPoint( DVector3( (**showItr).Apx_x,
                                          (**showItr).Apx_y,
                                          (**showItr).Apx_z ) );
        
        photon->setFitLayPointErr( DVector3( (**showItr).error_Apx_x,
                                             (**showItr).error_Apx_y,
                                             (**showItr).error_Apx_z ) );

        photon->setFitLaySlope( DVector3( (**showItr).Cx,
                                          (**showItr).Cy,
                                          (**showItr).Cz ) );
        
        photon->setFitLaySlopeErr( DVector3( (**showItr).error_Cx,
                                             (**showItr).error_Cy,
                                             (**showItr).error_Cz ) );
        
        _data.push_back( photon );
    }
    
    return NOERROR;
}
