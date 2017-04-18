/*
 *  GammaZToXYZ.cc
 *  GlueXTools
 *
 *  Modified GammaPToXYP.cc, replacing proton with heavy Z target
 *  Elton 4/14/2017
 *
 */

#include <iostream>
#include "TLorentzVector.h"

#include "AMPTOOLS_MCGEN/GammaZToXYZ.h"
#include "AMPTOOLS_MCGEN/TwoBodyDecayFactory.h"

#include "IUAmpTools/Kinematics.h"


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
    Double_t sigmagg = 1;         // take sigma (gg -> pi+ pi-) = 1 for now
    Double_t Z = 50;              // Z of Sn, target
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
        f = (coef/2)* (sigmagg/Wpipi)*Eg*Eg*Eg*Eg * (t0-t)* betapipi*betapipi * (FF*FF/(t*t))*conv*conv*conv*conv/(p1cm*p3cm*p1cm*p3cm);
    }
    else   {
        f = 0;
    }
    
    // cout << " t=" << t << " betastar=" << betastar << " gammastar=" << gammastar << " betapipicm=" << betapipicm << " t0=" << t0 << " thepipicm=" << thepipicm << " thepipi=" << thepipi << " f=" << f << endl;
    // cout << " t=" << t << " FF=" << FF << " f=" << f << endl;
    return	f;
}


// FORTRAN routines
extern "C"{
void cobrems_(float* Emax, float* Epeak, float* emitmr, float* radt, float* Dist, float* collDiam, int* doPolFluxfloat);
float dntdx_(float* x);
float dnidx_(float* x);
};

// Wrapper function for total
double dNtdx(double x)
{
        float xx = x;
        return (double)dntdx_(&xx);
}

double dNidx(double x)
{
        float xx = x;
        return (double)dnidx_(&xx);
}

GammaZToXYZ::GammaZToXYZ( float lowMassXY, float highMassXY, 
                          float massX, float massY, float beamMaxE, float beamPeakE, float beamLowE, float beamHighE,
                          ProductionMechanism::Type type ) : 

m_prodMech( ProductionMechanism::kZ, type, 6.0 ), // last arg is t dependence
// m_target( 0, 0, 0, 108.),    // use mass of Tin
m_target( 0, 0, 0, 208.),    // use mass of Pb since it is defined in particle tables.
m_childMass( 0 ) {

  m_childMass.push_back( massX );
  m_childMass.push_back( massY );
  
  m_prodMech.setMassRange( lowMassXY, highMassXY );
 
  // Initialize coherent brem table
  float Emax =  beamMaxE;
  float Epeak = beamPeakE;
  float Elow = beamLowE;
  float Ehigh = beamHighE;
  
  int doPolFlux=0;  // want total flux (1 for polarized flux)
  float emitmr=10.e-9; // electron beam emittance
  float radt=20.e-6; // radiator thickness in m
  float collDiam=0.0034; // meters
  float Dist = 76.0; // meters
  cobrems_(&Emax, &Epeak, &emitmr, &radt, &Dist, &collDiam, &doPolFlux);

  // Create histogram
  cobrem_vs_E = new TH1D("cobrem_vs_E", "Coherent Bremstrahlung vs. E_{#gamma}", 1000, Elow, Ehigh);
  
  // Fill histogram
  for(int i=1; i<=cobrem_vs_E->GetNbinsX(); i++){
	  double x = cobrem_vs_E->GetBinCenter(i)/Emax;
	  double y = 0;
	  if(Epeak<Elow) y = dNidx(x);
	  else y = dNtdx(x);
	  cobrem_vs_E->SetBinContent(i, y);
  }


    Int_t const npar = 4;
    Double_t Wpipi = 0.32;
    Double_t Eg=9;
    Double_t R0 = 5.358;   // Sn half-density radius, fm
    Double_t a0 = 0.550;   // Sn difuseness parameter, fm
    
    Double_t tmin=0;            // use -t (positive) as input
    Double_t tmax=0.05;
    Double_t xin[1];
    Double_t parin[npar];
    parin[0] = Wpipi;
    parin[1] = Eg;
    parin[2] = R0;
    parin[3] = a0;

  // Create tdist histogram
    Primakoff_tdist = new TH1D("Primakoff_tdist", "Primakoff t distribution", 1000, tmin,tmax);
  
  // Fill histogram
  for(int i=1; i<=Primakoff_tdist->GetNbinsX(); i++){
	  xin[0] = Primakoff_tdist->GetBinCenter(i);
	  double y = sigmat_func(xin,parin);
	  Primakoff_tdist->SetBinContent(i, y);
  }

}

Kinematics* 
GammaZToXYZ::generate(){

  double beamE = cobrem_vs_E->GetRandom();
  m_beam.SetPxPyPzE(0,0,beamE,beamE);

  double t = -Primakoff_tdist->GetRandom();       // histogram uses positive -t. 

  TLorentzVector resonance = m_prodMech.produceResonanceZ( m_beam , t);
  double genWeight = m_prodMech.getLastGeneratedWeight();
  
  vector< TLorentzVector > allPart;
  allPart.push_back( m_beam );
  allPart.push_back( m_beam + m_target - resonance );
  
  TwoBodyDecayFactory decay( resonance.M(), m_childMass );
  
  vector<TLorentzVector> fsPart = decay.generateDecay();
  
  for( vector<TLorentzVector>::iterator aPart = fsPart.begin();
      aPart != fsPart.end(); ++aPart ){
    
    aPart->Boost( resonance.BoostVector() );
    allPart.push_back( *aPart );
  }
 
  return new Kinematics( allPart, genWeight );
}

void
GammaZToXYZ::addResonance( float mass, float width, float bf ){
  
  m_prodMech.addResonance( mass, width, bf );
}

