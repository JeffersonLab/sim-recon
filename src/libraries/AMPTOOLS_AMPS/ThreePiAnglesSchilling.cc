
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/ThreePiAnglesSchilling.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

ThreePiAnglesSchilling::ThreePiAnglesSchilling( const vector< string >& args ) :
UserAmplitude< ThreePiAnglesSchilling >( args )
{
	assert( args.size() == 9 );
	
	rho000  = AmpParameter( args[0] );
	rho100  = AmpParameter( args[1] );
	rho1m10 = AmpParameter( args[2] );
  
	rho111  = AmpParameter( args[3] );
	rho001  = AmpParameter( args[4] );
	rho101  = AmpParameter( args[5] );
	rho1m11 = AmpParameter( args[6] );
	
	rho102  = AmpParameter( args[7] );
	rho1m12 = AmpParameter( args[8] );

	// need to register any free parameters so the framework knows about them
	registerParameter( rho000 );
	registerParameter( rho100 );
	registerParameter( rho1m10 );

	registerParameter( rho111 );
	registerParameter( rho001 );
	registerParameter( rho101 );
	registerParameter( rho1m11 );

	registerParameter( rho102 );
	registerParameter( rho1m12 );
}


complex< GDouble >
ThreePiAnglesSchilling::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
    TLorentzVector p3     ( pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0] );

	TLorentzVector resonance = p1 + p2 + p3;
	TLorentzRotation resonanceBoost( -resonance.BoostVector() );
	
	TLorentzVector beam_res = resonanceBoost * beam;
	TLorentzVector recoil_res = resonanceBoost * recoil;
	TLorentzVector p1_res = resonanceBoost * p1;
    TLorentzVector p2_res = resonanceBoost * p2;
	
	TVector3 z = -recoil_res.Vect().Unit();
	TVector3 y = beam_res.Vect().Cross(z).Unit();
	TVector3 x = y.Cross(z).Unit();

    // For the three pion case, we define the angles relative to the normal of the decay plane
    TVector3 norm = p1_res.Vect().Cross(p2_res.Vect()).Unit();
	
	TVector3 angles( norm.Dot(x), norm.Dot(y), norm.Dot(z) );
	
	GDouble cosTheta = angles.CosTheta();
	GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
	GDouble sin2Theta = sin(2.*angles.Theta());
	
	GDouble phi = angles.Phi();
	GDouble Phi = recoil.Vect().Phi();
	
	//GDouble psi = phi - Phi;
	//if(psi < -1*PI) psi += 2*PI;
	//if(psi > PI) psi -= 2*PI;
	
	// vector meson production from K. Schilling et. al.
	GDouble Pgamma = 0.5;
	
	GDouble W = 0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*cosTheta*cosTheta - sqrt(2.)*rho100*sin2Theta*cos(phi) - rho1m10*sinSqTheta*cos(2.*phi);
	
	W -= Pgamma*cos(2.*Phi) * (rho111*sinSqTheta + rho001*cosTheta*cosTheta - sqrt(2.)*rho101*sin2Theta*cos(phi) - rho1m11*sinSqTheta*cos(2.*phi));
	
	W -= Pgamma*sin(2.*Phi) * (sqrt(2.)*rho102*sin2Theta*sin(phi) + rho1m12*sinSqTheta*sin(2.*phi));
	
	W *= 3./(4.*PI);

	return complex< GDouble > ( sqrt(fabs(W)) );
}

