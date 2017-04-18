
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiW_primakoff.h"

// Class modeled after BreitWigner amplitude function provided for examples with AmpTools.
// Dependence of swave 2pi cross section on W (mass of 2pi system)
// Elton 4/17/2017

TwoPiW_primakoff::TwoPiW_primakoff( const vector< string >& args ) :
UserAmplitude< TwoPiW_primakoff >( args )
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
TwoPiW_primakoff::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector P1, P2, Ptot, Ptemp;
  
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
  
  GDouble mass  = Ptot.M();
  GDouble mass1 = P1.M();
  GDouble mass2 = P2.M();

  cout << "calcAmplitude: mass=" << mass << " mass1=" << mass1 << " mass2=" << mass2 << endl;
  
  complex<GDouble> bwtop( mass1, 0.0 );
  
  complex<GDouble> bwbottom( mass2, 0.0 );
  
  return( bwtop / bwbottom );
}

void
TwoPiW_primakoff::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


