
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Pi0Regge.h"

// FORTRAN routines
extern "C"{
void test_(float* Eg, float* t, float* DSG);
void diffcross_(float* Eg, float* t, float* DSG);
};

// Wrapper function for cross section
GDouble DiffCross(double x, double y)
{
        float xx = x;
	float yy = y;
	float zz = 0;
        diffcross_(&xx, &yy, &zz);

	return (double)zz;
}


Pi0Regge::Pi0Regge( const vector< string >& args ) :
UserAmplitude< Pi0Regge >( args )
{

	/*
	assert( args.size() == 9 );
	
	rho000  = AmpParameter( args[0] );

	// need to register any free parameters so the framework knows about them
	registerParameter( rho000 );
	*/
}


complex< GDouble >
Pi0Regge::calcAmplitude( GDouble** pKin ) const {
  
  TLorentzVector target  ( 0., 0., 0., 0.938);
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 

  TLorentzVector cm = recoil + p1;
  TLorentzRotation cmBoost( -cm.BoostVector() );

  TLorentzVector target_cm = cmBoost * target;
  TLorentzVector recoil_cm = cmBoost * recoil;

  // some variables for photon polarization dependence (not currently implemented)
  //TLorentzVector p1_cm = cmBoost * p1;
  //GDouble phi = p1_cm.Phi();
  //if(phi < -1*PI) phi += 2*PI;
  //if(phi > PI) phi -= 2*PI;
  //GDouble cosTheta = p1_cm.CosTheta();
  //GDouble Pgamma = 0.5;
  
  // factors needed to calculate amplitude in fortran code
  float Eg = beam.E();
  float t = (target_cm - recoil_cm).M2();

  GDouble W = DiffCross(Eg, t);

  return complex< GDouble > ( sqrt(W) );
}

