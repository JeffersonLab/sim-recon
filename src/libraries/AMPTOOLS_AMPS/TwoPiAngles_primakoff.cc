
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
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
	TLorentzVector resonance = p1 + p2;
  
	
	// production plane is defined by the pi+ (neglect recoil)
	TVector3 zlab(0.,0.,1.0);     // z axis in lab
	TVector3 y = (p1.Vect().Cross(zlab)).Unit();    // perpendicular to decay plane. ensure that y is perpendicular to z

        TVector3 eps(cos(phipol), sin(phipol), 0.0); // beam polarization vector in lab
	TVector3 eps_perp = eps.Cross(zlab).Unit();         // perpendicular to plane defined by eps
        GDouble Phi_pip = atan2(y.Dot(eps),y.Dot(eps_perp));  // use this calculation to preserve sign of angle
	// GDouble Phi2 = acos(y.Dot(eps_perp));
	/*cout << endl << "phipol=" << phipol << " p1.Phi=" << p1.Phi() << endl;
	cout << " p1 Angles="; p1.Vect().Print();
	cout << "y= "; y.Print();
	cout << "eps_perp= "; eps_perp.Print();
	cout << "eps= "; eps.Print();*/
	// Phi_pip = Phi_pip > 0? Phi_pip : Phi_pip + 3.14159;                     // make angle between eps and decay plane a positive number. This is psi of vector-meson production in forward kinematics.

        // Get cosTheta and phi using helicity frame, but should be irrelevant for the s-wave Primakoff. Still use y from plane of 2 pions
	// May be useful when rho production is added as background.
	TLorentzRotation resonanceBoost( -resonance.BoostVector() );
	
	TLorentzVector beam_res = resonanceBoost * beam;
	TLorentzVector recoil_res = resonanceBoost * recoil;
	TLorentzVector p1_res = resonanceBoost * p1;

        // choose helicity frame: z-axis opposite recoil target in rho rest frame. Note that for Primakoff recoil is never measured.
        y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();   // redefine y for self-consistency
        TVector3 z = -1. * recoil_res.Vect().Unit();
        TVector3 x = y.Cross(z).Unit();
        TVector3 angles( (p1_res.Vect()).Dot(x),
                         (p1_res.Vect()).Dot(y),
                         (p1_res.Vect()).Dot(z) );

        GDouble CosTheta = angles.CosTheta();
        // GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
        // GDouble sin2Theta = sin(2.*angles.Theta());

	// calculate Phi based on recoil momenta (of course not measureable for Primakoff events)

        // GDouble Phi_prod = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
	// Phi_prod = Phi_prod > 0? Phi_prod : Phi_prod + 3.14159;

        GDouble phi = angles.Phi();
        GDouble psi = Phi_pip;               // in the limit of forward scattering (Primakoff), Phi_pip is the angle between pip and the polarization
	// GDouble Phi = Phi_prod;              // retain angle between hadronic plane and polarization, although random for Primakoff
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
	  prefactor = 0.5*sqrt(1-polFrac)*(cos(Phi_pip)-i*sin(Phi_pip) - eta_c*(cos(Phi_pip)+i*sin(Phi_pip)) );
	  Mrho = m_rho;
	  break;
        case 1:
	  prefactor = 0.5*sqrt(1+polFrac)*(cos(Phi_pip)-i*sin(Phi_pip) + eta_c*(cos(Phi_pip)+i*sin(Phi_pip)) );
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

