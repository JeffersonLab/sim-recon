
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_primakoff.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

TwoPiAngles_primakoff::TwoPiAngles_primakoff( const vector< string >& args ) :
UserAmplitude< TwoPiAngles_primakoff >( args )
{
	assert( args.size() == 5 );
	
	phipol  = atof(args[0].c_str() )*3.14159/180.; // azimuthal angle of the photon polarization vector in the lab. Convert to radians.
	polFrac  = AmpParameter( args[1] ); // fraction of polarization (0-1)
	m_rho = atoi( args[2].c_str() );  // Jz component of rho 
	PhaseFactor  = AmpParameter( args[3] );  // prefix factor to amplitudes in computation
	flat = atoi( args[4].c_str() );  // flat=1 uniform angles, flat=0 use YLMs 

	assert( ( phipol >= 0.) && (phipol <= 2*3.14159));
	assert( ( polFrac >= 0 ) && ( polFrac <= 1 ) );
        assert( ( m_rho == 0 ) );
        assert( ( PhaseFactor == 0 ) || ( PhaseFactor == 1 ));
	assert( (flat == 0) || (flat == 1) );

	// need to register any free parameters so the framework knows about them
	registerParameter( polFrac );
}


complex< GDouble >
TwoPiAngles_primakoff::calcAmplitude( GDouble** pKin ) const {

  // for Primakoff, all calculations are in the lab frame. Keep recoil but remember that it cannot be measured by detector.
  
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
	TLorentzVector p1     ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p2     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector recoil ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
	TLorentzVector resonance = p1 + p2;

        TVector3 eps(cos(phipol), sin(phipol), 0.0); // beam polarization vector in lab

	TLorentzRotation resonanceBoost( -resonance.BoostVector() );
	
	TLorentzVector beam_res = resonanceBoost * beam;
	TLorentzVector recoil_res = resonanceBoost * recoil;
	TLorentzVector p1_res = resonanceBoost * p1;

        // choose helicity frame: z-axis opposite recoil target in rho rest frame. Note that for Primakoff recoil is defined as missing P4
        TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();   
        TVector3 z = -1. * recoil_res.Vect().Unit();
        TVector3 x = y.Cross(z).Unit();
        TVector3 angles( (p1_res.Vect()).Dot(x),
                         (p1_res.Vect()).Dot(y),
                         (p1_res.Vect()).Dot(z) );

        GDouble CosTheta = angles.CosTheta();
        GDouble phi = angles.Phi();
        // GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
        // GDouble sin2Theta = sin(2.*angles.Theta());

        GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

        GDouble psi = Phi - phi;               // define angle difference 
        if(psi < -1*PI) psi += 2*PI;
        if  (psi > PI) psi -= 2*PI;

	/*cout << " recoil_res Angles="; recoil_res.Vect().Print();
	cout << " p1_res Angles="; p1_res.Vect().Print();
	cout << "Phi_pip= " << Phi_pip << endl;
	cout << "Phi= " << Phi << endl;
	cout << "Phi_prod= " << Phi_prod << endl;
	cout << "phi= " << phi << endl;
	cout << " psi=" << psi << endl;*/
     
	complex< GDouble > i( 0, 1 );
	complex< GDouble > prefactor( 0, 0 );
	complex< GDouble > Amp( 0, 0 );
	complex< GDouble> eta_c(1,0);
	Int_t Mrho=0;

	switch (PhaseFactor) {
        case 0:
	  prefactor = 0.5*sqrt(1-polFrac)*(cos(Phi)-i*sin(Phi) - eta_c*(cos(Phi)+i*sin(Phi)) );
	  Mrho = m_rho;
	  break;
        case 1:
	  prefactor = 0.5*sqrt(1+polFrac)*(cos(Phi)-i*sin(Phi) + eta_c*(cos(Phi)+i*sin(Phi)) );
	  Mrho = m_rho;
	  break;
	}
	
	if (flat == 1) {
	  Amp = 1;
	}
	else {
	  Amp =  prefactor * Y( 0, Mrho, CosTheta, phi);
	}

	// cout << " m_rho=" << m_rho << " CosTheta=" << CosTheta << " phi=" << phi << " prefactor=" << prefactor << " Amp=" << Amp << endl;

	return Amp;
}

