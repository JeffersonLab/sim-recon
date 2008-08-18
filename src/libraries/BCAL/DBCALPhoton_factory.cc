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

  //scaling parameter set for nCell < 5
   
    m_scaleZ_p0LT =  7.038E-01;
    m_scaleZ_p1LT =  1.251E-03;
    m_scaleZ_p2LT = -2.539E-06;
     
    m_nonlinZ_p0LT =  1.181E-01;
    m_nonlinZ_p1LT = -9.330E-04;
    m_nonlinZ_p2LT =  4.334E-06;    
    m_nonlinZ_p3LT = -5.486E-09;

   //scaling parameter set for nCell >= 5
   
    m_scaleZ_p0GE =  8.483E-01;
    m_scaleZ_p1GE =  5.840E-04;
    m_scaleZ_p2GE = -1.690E-06;
     
    m_nonlinZ_p0GE =  2.489E-02;
    m_nonlinZ_p1GE =  4.504E-04;
    m_nonlinZ_p2GE = -1.980E-06;    
    m_nonlinZ_p3GE =  2.428E-09;

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
    loop->Get( showerVect );
    
    for( vector< const DBCALShower* >::iterator showItr = showerVect.begin();
         showItr != showerVect.end();
         ++showItr ){
                
        double xSh = (**showItr).x;
        double ySh = (**showItr).y;
        double zSh = (**showItr).z;
        int nCell = (**showItr).N_cell;        

        // get z where shower enters BCAL (this corresponds to generated z
        // in tuning MC)
        double zEntry = zSh - ( ( zSh - m_zTarget ) * 
                                ( 1 - m_bcalIR / ( sqrt( xSh * xSh + ySh * ySh ) ) ) ); 
        
        // calibrate energy:
        // Energy calibration has a z dependence -- the
        // calibration comes from fitting E_rec / E_gen to scale * E_gen^nonlin
        // for slices of z.  These fit parameters (scale and nonlin) are then plotted 
        // as a function of z and fit.
        
	if( nCell < 5 ) { 
            nonlin = m_nonlinZ_p0LT +
            m_nonlinZ_p1LT * zEntry +
            m_nonlinZ_p2LT * zEntry * zEntry +
	    m_nonlinZ_p3LT * zEntry * zEntry * zEntry;}
	else if( nCell >= 5 ) {
            nonlin = m_nonlinZ_p0GE +
            m_nonlinZ_p1GE * zEntry +
            m_nonlinZ_p2GE * zEntry * zEntry +
	    m_nonlinZ_p3GE * zEntry * zEntry * zEntry;}
	//	nonlin = 0.0; // fixed value for debug

        // for the scale, we extend the scale factor at the target as a constant behind
        // the target
        if( zEntry < m_zTarget ) zEntry = m_zTarget;
    
        if( nCell < 5 ) {
            scale = m_scaleZ_p0LT +
            m_scaleZ_p1LT * zEntry +
	    m_scaleZ_p2LT * zEntry * zEntry;}
	else if ( nCell >= 5 ) {
            scale = m_scaleZ_p0GE +
            m_scaleZ_p1GE * zEntry +
	    m_scaleZ_p2GE * zEntry * zEntry;}
	//     scale = 1.0; // fixed value for debug
        
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
