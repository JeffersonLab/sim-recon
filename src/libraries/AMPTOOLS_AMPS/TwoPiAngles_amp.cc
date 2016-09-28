
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_amp.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

TwoPiAngles_amp::TwoPiAngles_amp( const vector< string >& args ) :
UserAmplitude< TwoPiAngles_amp >( args )
{
	assert( args.size() == 4 );
	
	phipol  = atof(args[0].c_str() )*3.14159/180.; // azimuthal angle of the photon polarization vector in the lab.
	polFrac  = AmpParameter( args[1] ); // fraction of polarization (0-1)
	m_rho = atoi( args[2].c_str() );  // Jz component of rho 
	PhaseFactor  = AmpParameter( args[3] );  // prefix factor to amplitudes in computation ( 0=1/1=exp(2iPhi)/2=-exp(2iPhi) )



	assert( ( phipol >= 0.) && (phipol <= 2*3.14159));
	assert( ( polFrac >= 0 ) && ( polFrac <= 1 ) );
        assert( ( m_rho == 1 ) || ( m_rho == 0 ) || ( m_rho == -1 ));
        assert( ( PhaseFactor == 0 ) || ( PhaseFactor == 1 ) || ( PhaseFactor == 2 ) || ( PhaseFactor == 3 ));



	// need to register any free parameters so the framework knows about them
	registerParameter( polFrac );
}


complex< GDouble >
TwoPiAngles_amp::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
	
	TLorentzVector resonance = p1 + p2;
	TLorentzRotation resonanceBoost( -resonance.BoostVector() );
	
	TLorentzVector beam_res = resonanceBoost * beam;
	TLorentzVector recoil_res = resonanceBoost * recoil;
	TLorentzVector p1_res = resonanceBoost * p1;
	
	// normal to the production plane
        TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

        // choose helicity frame: z-axis opposite recoil proton in rho rest frame
        TVector3 z = -1. * recoil_res.Vect().Unit();
        TVector3 x = y.Cross(z).Unit();
        TVector3 angles( (p1_res.Vect()).Dot(x),
                         (p1_res.Vect()).Dot(y),
                         (p1_res.Vect()).Dot(z) );

        GDouble cosTheta = angles.CosTheta();
        // GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
        // GDouble sin2Theta = sin(2.*angles.Theta());
        GDouble phi = angles.Phi();

	// TVector3 zlab(0.,0.,1.0);     // z axis in lab
        TVector3 eps(cos(phipol), sin(phipol), 0.0); // beam polarization vector in lab
	// TVector3 eps_perp = zlab.Cross(eps);         // perpendicular to plane defined by eps
	// GDouble Phi_test = asin((eps_perp.Cross(y)).Mag());        // compute angle between planes. 
        GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
	Phi = Phi > 0? Phi : Phi + 3.14159;

	// cout << "Phi_test=" << Phi_test << " Phi=" << Phi << " Sum=" << Phi_test+Phi << " Diff=" << Phi_test-Phi << " PhaseFactor=" << PhaseFactor << endl;
     
	complex< GDouble > i( 0, 1 );
	complex< GDouble > prefactor( 0, 0 );
	complex< GDouble > Amp( 0, 0 );

	switch (PhaseFactor) {
        case 0:
	  prefactor = 0.5*sqrt(1-polFrac);
	  break;
        case 1:
	  prefactor = 0.5*sqrt(1+polFrac);
	  break;
        case 2:
	  prefactor = 0.5*sqrt(1-polFrac)*(cos(2*Phi) + i*sin(2*Phi));
	  break;
        case 3:
	  prefactor = -0.5*sqrt(1+polFrac)*(cos(2*Phi) - i*sin(2*Phi));
          break;
	}
	
	Amp =  prefactor * Y( 1, m_rho, cosTheta, phi);

	// cout << " m_rho=" << m_rho << " cosTheta=" << cosTheta << " phi=" << phi << " prefactor=" << prefactor << " Amp=" << Amp << endl;

	return Amp;
}

