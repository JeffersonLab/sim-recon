
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiWt_primakoff.h"

// Class modeled after BreitWigner amplitude function provided for examples with AmpTools.
// Dependence of swave 2pi cross section on W (mass of 2pi system)
// Elton 4/17/2017


Double_t sigma_ggpipi_func (Double_t *x, Double_t *par){
    
    // Parameterization from Rory for cross section for gamma gamma -> pi+ pi-
    // Returns cross section in units of nb/sr
    
    // constants
    // Double_t const PI = 3.14159;
    Double_t MPI =0.139570;
    Double_t W0 = 0.3;
    
    Double_t expon = par[0];
    // Double_t par2 = par[1];
    Double_t Wpipi = x[0] ;
    Double_t f;
    
    if (Wpipi < 2*MPI) {
        f = 0;
    }
    else if (Wpipi < W0) {
        f = 300./(0.6*4.*PI)*(Wpipi-2.*MPI)/(W0-2.*MPI);           // linear rise, isotropic, CB data only 60% coverage
    }
    else {
        f = 300./(0.6*4.*PI)*pow(W0/Wpipi,expon);	                 // power fall off, isotropic
    }
    
    return	f;
}


Double_t ff_func (Double_t *x, Double_t *par){
    
    // return the nuclear form factor accourding to 2 parameter Fermi distribution
    // See Journall of Research of the National Bureau of Standards - B Mathenatics and Mathematical Physics
    // Vol. 70B, No. 1, Jan-Mar 1966. "The Form Factor of the Fermi Model Spatial Distribution," by Maximon and Schrack
    //
    // Function is a function of q, the three-momentum transfer to the nucleus.
    // Input argument is t
    // Note that q is the 3-vector momentum, but for low -t, q ~ sqrt(-t).
    
    // constants
    // Double_t alpha = 1/137.;
    Double_t pi = 3.14159;
    Double_t hbarc = 0.19733;                  // GeV*fm
    Double_t q = sqrt(x[0])/hbarc;          // take q to be  in fm^-1. Input variable is positive (-t)
    
    Double_t R0  = par[0];  				 // half-density radius
    Double_t a0 = par[1];                    // skin or diffuseness parameter
    
    Double_t rho0;
    Double_t sum=0;
    Int_t jmax=4;
    for (Int_t j=1;j<jmax;j++) {                    // truncate after 3 terms, first term dominates.
        Double_t sum1 =pow(-1.,j-1)*exp(-j*R0/a0)/(j*j*j);
        sum += sum1;
        // cout << "jmax=" << jmax << " j=" << j << " R0=" << R0 << " a0=" << a0 << " sum1=" << sum1 << " sum=" << sum << endl;
    }
    
    rho0 = (4*pi*R0/3)*( pi*a0*pi*a0  + R0*R0 + 8*pi*a0*a0*a0*sum);
    rho0 = 1/rho0;
    
    Double_t f = 0;
    
    f = (4*pi*pi*rho0*a0*a0*a0)/(q*q*a0*a0*sinh(pi*q*a0)*sinh(pi*q*a0))
    	* (pi*q*a0 * cosh(pi*q*a0) * sin(q*R0) - q*R0 *cos(q*R0) * sinh(pi*q*a0) )
    	+ 8*pi*rho0*a0*a0*a0*sum;
    
    // cout << " q=" << q << " f=" << f << endl;
    return f;
    
}


Double_t sigmat_func (Double_t *x, Double_t *par){
    
    // return the cross section for Primakoff production of pi+pi-. CPP proposal PR12-13-008 Eq 8.
    // independent variable is momentum transfer -t, in GeV2;
    
    // constants
    Double_t alpha = 1/137.;
    Double_t pi = 3.14159;
    Double_t betapipi = 0.999;      // beta=0.999 (W=0.3), beta=0.997 (W=0.4), beta= 0.986 (W=1.0 GeV)
    // Double_t sigmagg = 1;         // take sigma (gg -> pi+ pi-) = 1 for now
    // Double_t Z = 50;              // Z of Sn, target
    Double_t Z = 82;              // Z of Pb, target
    Double_t coef = 4*alpha*Z*Z/(pi);
    
    Double_t Wpipi = par[0];
    Double_t Eg = par[1];
    Double_t t = -x[0] ;     // theta of pipi system in the lab. Input Variable is positive (i.e. -t)
    
    Double_t xin[2];
    Double_t parin[2];
    
    xin[0] = -t;                       // input variable to ff is positive (-t)
    parin[0] = par[2];  	       // half-density radius, fm
    parin[1] = par[3];                 // diffuseness paramter, fm
    
    Double_t FF = ff_func(xin,parin);
    
    // include other masses here
    
    Double_t m1 = 0;      // mass of photon, incident beam
    // Double_t m2 = 108.;    // mass of 116Sn, target
    Double_t m2 = 208.;    // use Pb mass because it is in the particle list
    Double_t m3 = Wpipi;   // mass of 2pi system, scattered system
    Double_t m4 = m2;     // recoil target
    
    
    Double_t f = 0;
    
    Double_t s = m2*m2 + 2*Eg*m2;
    Double_t sqrts = s > 0? sqrt(s) : 0;
    if (s < 0) {
        cout << "*** sigma_func: s =" << s << " < 0!" << endl;
        return f;
    }
   
    Double_t E1cm = (s+m1*m1-m2*m2)/(2*sqrts);
    // Double_t E2cm = (s+m2*m2-m1*m1)/(2*sqrts);
    Double_t E3cm = (s+m3*m3-m4*m4)/(2*sqrts);
    // Double_t E4cm = (s+m4*m4-m3*m3)/(2*sqrts);
    
    Double_t p1cm = E1cm*E1cm - m1*m1? sqrt(E1cm*E1cm - m1*m1) : 0;
    // Double_t p2cm = E2cm*E2cm - m2*m2? sqrt(E2cm*E2cm - m2*m2) : 0;
    Double_t p3cm = E3cm*E3cm - m3*m3? sqrt(E3cm*E3cm - m3*m3) : 0;
    // Double_t p4cm = E4cm*E4cm - m4*m4? sqrt(E4cm*E4cm - m4*m4) : 0;
    
    Double_t arg = (m1*m1-m3*m3-m2*m2+m4*m4)/(2*sqrts);
    Double_t t0 = arg*arg - (p1cm - p3cm)*(p1cm - p3cm);
    // Double_t t1 = arg*arg - (p1cm + p3cm)*(p1cm + p3cm);
    
    Double_t betastar = Eg/(Eg + m2);
    Double_t gammastar = (Eg + m2)/sqrts;
    Double_t betapipicm = p3cm/E3cm;
    
    // Double_t thepipicm = t0 -t > 0? sqrt((t0 -t)/(p1cm*p3cm)) : 0;
    
    Double_t conv = 1./(gammastar*(1 + betastar/betapipicm));
    
    if (-t > -t0) {
        f = (coef/2)* Eg*Eg*Eg*Eg * (t0-t)* betapipi*betapipi * (FF*FF/(t*t))*conv*conv*conv*conv/(p1cm*p3cm*p1cm*p3cm);
    }
    else   {
        f = 0;
    }
    
    // cout << " t=" << t << " betastar=" << betastar << " gammastar=" << gammastar << " betapipicm=" << betapipicm << " t0=" << t0 << " thepipicm=" << thepipicm << " thepipi=" << thepipi << " f=" << f << endl;
    // cout << " t=" << t << " FF=" << FF << " f=" << f << endl;
    return	f;
}


TwoPiWt_primakoff::TwoPiWt_primakoff( const vector< string >& args ) :
UserAmplitude< TwoPiWt_primakoff >( args )
{
  
  assert( args.size() == 4 );
	m_par1 = AmpParameter( args[0] );
	m_par2 = AmpParameter( args[1] );
	m_daughters = pair< string, string >( args[2], args[3] );
  
  // need to register any free parameters so the framework knows about them
  registerParameter( m_par1 );
  registerParameter( m_par2 );
  
  // make sure the input variables look reasonable
  // assert( ( m_orbitL >= 0 ) && ( m_orbitL <= 4 ) );
}

complex< GDouble >
TwoPiWt_primakoff::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector P1, P2, Ptot, Ptemp, Precoil;
  
  for( unsigned int i = 0; i < m_daughters.first.size(); ++i ){
    
    string num; num += m_daughters.first[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    P1 += Ptemp;
    Ptot += Ptemp;
  }
  
  for( unsigned int i = 0; i < m_daughters.second.size(); ++i ){
    
    string num; num += m_daughters.second[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    P2 += Ptemp;
    Ptot += Ptemp;
  }
  
  GDouble Wpipi  = Ptot.M();
  // GDouble mass1 = P1.M();
  // GDouble mass2 = P2.M();

  // get momentum transfer
  Precoil.SetPxPyPzE (pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);
  GDouble Et = Precoil.E();
  GDouble Mt = Precoil.M();
  GDouble t = -2*Precoil.M()*(Et - Mt);      

  // cout << "Precoil.M()=" << Precoil.M() << " T=" << Precoil.E() - Precoil.M() << " t=" << t << endl; Precoil.Print();

  // call sigma (gamma gamma -> pi pi) cross section

  
    Int_t const npar = 4;
    Double_t xin[1];
    xin[0] = Wpipi;                // W, 2pi mass
    Double_t Eg = pKin[0][0];          // incident photon energy
    Double_t parin[npar];
    parin[0] = 1.29;              // parameter 1: exponent
    parin[1] = 0.;                // parameter 2: par2 (spare)
    // Double_t Wmin=0.2 ;
    // Double_t Wmax=0.8;

    GDouble sig_ggpipi = sigma_ggpipi_func(xin,parin);

    parin[0] = Wpipi;
    parin[1] = Eg;
    // Double_t R0 = 6.62;   // Pb half-density radius, fm
    // Double_t a0 = 0.546;   // Pb difuseness parameter, fm
    // Double_t R0  = 5.358;   // Sn half-density radius, fm
    // Double_t a0 = 0.550;   // Sn difuseness parameter, fm
    parin[2] = 6.62;   
    parin[3] = 0.546;  
    xin[0] = -t;             // input positive value of t
    
    GDouble sigmat = sigmat_func (xin,parin);

    // cout << "calcAmplitude: 2pi mass=" << Wpipi << " Eg=" << Eg << " t=" << t << " sig_ggpipi=" << sig_ggpipi << " sigmat=" << sigmat << endl;
  
    complex<GDouble> Csig( sqrt(sigmat*sig_ggpipi/Wpipi/exp(6.0*t)), 0.0 );    // Return complex double, sqrt (cross section). Divide out generated exponential. 
  
  return( Csig  );
}

void
TwoPiWt_primakoff::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


