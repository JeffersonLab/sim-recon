
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Pi0Regge.h"

Pi0Regge::Pi0Regge( const vector< string >& args ) :
UserAmplitude< Pi0Regge >( args )
{
	assert( args.size() == 2 );
	Pgamma = atof( args[0].c_str() );
	// Polarization plane angle (PARA = 0 and PERP = PI/2)
	PolPlane = atof( args[1].c_str() );
}


complex< GDouble >
Pi0Regge::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector target  ( 0., 0., 0., 0.938);
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	
	TLorentzVector cm = recoil + p1;
	TLorentzRotation cmBoost( -cm.BoostVector() );
	
	TLorentzVector beam_cm = cmBoost * beam;
	TLorentzVector target_cm = cmBoost * target;
	TLorentzVector recoil_cm = cmBoost * recoil;
	
	// phi dependence needed for polarized distribution
	TLorentzVector p1_cm = cmBoost * p1;
	GDouble phi = p1_cm.Phi() + PolPlane*TMath::Pi()/180. + TMath::Pi()/2.;
	GDouble cos2Phi = cos(2.*phi);

	// factors needed to calculate amplitude in c++ code
	GDouble Ecom = cm.M();
	GDouble theta = p1_cm.Theta();	

	// amplitude coded in c++ (include calculation of beam asymmetry)
	GDouble BeamSigma = 0.;
	GDouble W = Pi0PhotCS_S(Ecom, theta, BeamSigma);
	W *= (1 - Pgamma * BeamSigma * cos2Phi);

	return complex< GDouble > ( sqrt( fabs(W) ) );
}

