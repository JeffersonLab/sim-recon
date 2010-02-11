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
//#include <iostream>
//using namespace std;

DBCALPhoton_factory::DBCALPhoton_factory()
{

  /*    
    m_scaleZ_p0GE =  8.25618e-01;
    m_scaleZ_p1GE =  9.22446e-02;
    m_scaleZ_p2GE =  1.27379e02;
    m_scaleZ_p3GE = -1.21171e02;
     
    m_nonlinZ_p0GE =  5.0833e-02;
    m_nonlinZ_p1GE = -2.66106e-02;
    m_nonlinZ_p2GE =  5.22213e01;    
    m_nonlinZ_p3GE =  6.88344e01;
    
    m_linZ_p0GE = -5.26475e-03;
    m_linZ_p1GE = -2.47419e-02;
    m_linZ_p2GE =  4.19082e01;    
    m_linZ_p3GE =  6.69810e01;
    
   
   //scaling parameter set for Z>370 (end of module)

    m_scaleZ_p0 =  0.8284;
    m_scaleZ_p1 =  3.3;
    m_scaleZ_p2 =  422.5;
    m_scaleZ_p3 =  12.04;
    m_scaleZ_p4 = 0.0;
     
    m_nonlinZ_p0 =  0.05136;
    m_nonlinZ_p1 = 1000.0;
    m_nonlinZ_p2 =  453.6;    
    m_nonlinZ_p3 = 17.21;

    m_linZ_p0 = -7.166e-03;
    m_linZ_p1 = -1000.0;
    m_linZ_p2 = 482.0;
    m_linZ_p3 = 24.19;
  */


 // parameters for correcting with dark noise (different correction function)

    m_scaleZ_p0GE =  0.94795;
    m_scaleZ_p1GE =  -1.747;
    m_scaleZ_p2GE =  1000;
    m_scaleZ_p3GE = 251.5;
     
    m_nonlinZ_p0GE =  0.02597;
    m_nonlinZ_p1GE = -0.0347;
    m_nonlinZ_p2GE =  54.89;    
    m_nonlinZ_p3GE =  55.53;
    
    m_linZ_p0GE = -3.0533e-03;
    m_linZ_p1GE = -0.01671;
    m_linZ_p2GE =  14.363;    
    m_linZ_p3GE =  69.542;
    
   
   //scaling parameter set for Z>370 (end of module)

    m_scaleZ_p0 =  0.8776;
    m_scaleZ_p1 =  -10.0;
    m_scaleZ_p2 =  428.2;
    m_scaleZ_p3 =  12.22;
    m_scaleZ_p4 = 0.0;
     
    m_nonlinZ_p0 = 0.02604;
    m_nonlinZ_p1 = -10.1;
    m_nonlinZ_p2 =  433;    
    m_nonlinZ_p3 = 13.38;

    m_linZ_p0 = -3.3587e-03;
    m_linZ_p1 = -10.0;
    m_linZ_p2 = 456.7;
    m_linZ_p3 = 17.69;

}

//------------------
// brun
//------------------
jerror_t DBCALPhoton_factory::brun(JEventLoop *loop, int runnumber)
{
    
    vector<const DBCALGeometry*> bcalGeomVect;
    loop->Get( bcalGeomVect );
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
	//     int nCell = (**showItr).N_cell;        

        // get z where shower enters BCAL (this corresponds to generated z
        // in tuning MC)
        double zEntry = zSh - ( ( zSh - m_zTarget ) * 
                                ( 1 - m_bcalIR / ( sqrt( xSh * xSh + ySh * ySh ) ) ) ); 
        
        // calibrate energy:
        // Energy calibration has a z dependence -- the
        // calibration comes from fitting E_rec / E_gen to scale * E_gen^nonlin
        // for slices of z.  These fit parameters (scale and nonlin) are then plotted 
        // as a function of z and fit.
    
  if( zEntry < 370.0 ) {
 scale = (m_scaleZ_p0GE  + m_scaleZ_p1GE *(exp( -0.5 *(zEntry - m_scaleZ_p2GE )* (zEntry - m_scaleZ_p2GE ) / (m_scaleZ_p3GE * m_scaleZ_p3GE)   ) ) );

  nonlin =( m_nonlinZ_p0GE  + m_nonlinZ_p1GE *(exp( -0.5 *(zEntry - m_nonlinZ_p2GE )* (zEntry - m_nonlinZ_p2GE ) / (m_nonlinZ_p3GE * m_nonlinZ_p3GE)   ) ) ) ;


  lin = ( m_linZ_p0GE  + m_linZ_p1GE *(exp( -0.5 *(zEntry - m_linZ_p2GE )* (zEntry - m_linZ_p2GE ) / (m_linZ_p3GE * m_linZ_p3GE)   ) ) ) ;

	//	nonlin = 0.0; // fixed value for debug
	//	   lin = 0.0; // fixed value for debug
   	//       scale = 1.0; // fixed value for debug
  }
  
        if( zEntry >= 370.0 ) {
            scale = m_scaleZ_p0 +  m_scaleZ_p1 *(exp( -0.5 *(zEntry - m_scaleZ_p2 )* (zEntry - m_scaleZ_p2 ) / (m_scaleZ_p3 * m_scaleZ_p3)   ) ) ;

  nonlin = m_nonlinZ_p0  + m_nonlinZ_p1 *(exp( -0.5 *(zEntry - m_nonlinZ_p2 )* (zEntry - m_nonlinZ_p2 ) / (m_nonlinZ_p3 * m_nonlinZ_p3)   ) )  ;


  lin = m_linZ_p0 + m_linZ_p1 *(exp( -0.5 *(zEntry - m_linZ_p2 )* (zEntry - m_linZ_p2 ) / (m_linZ_p3 * m_linZ_p3)   ) )  ;


 //  cout << scale << ' ' << nonlin << ' ' << lin << endl;    
	}


	// if( zEntry < m_zTarget ) zEntry = m_zTarget;

        //end of BCAL calibration        


        
        // now turn E_rec into E_gen -->> E_gen = ( E_rec / scale ) ^ ( 1 / ( 1 + nonlin ) )
        double energy = pow( ((**showItr).E - lin ) / scale, 1 / ( 1 + nonlin ) );
        
        
        double pScale = energy / 
            sqrt( xSh * xSh + ySh * ySh + ( zSh - m_zTarget ) * ( zSh - m_zTarget ) );
        
        DBCALPhoton* photon = new DBCALPhoton( (**showItr).id );
        
        photon->setLorentzMomentum( DLorentzVector( xSh * pScale, 
                                                    ySh * pScale, 
                                                    ( zSh - m_zTarget ) * pScale, 
                                                    energy ) );
        
        photon->setShowerPosition( DVector3( xSh, ySh, zSh ) );

        photon->setShowerTime( (**showItr).t );
        
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
