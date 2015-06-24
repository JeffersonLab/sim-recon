/*
  Generates N-body phase space decays using the Raubold Lynch
  method (F.James CERN 68-15 (1968).
  Borrows liberally from ROOT class TGenPhaseSpace as well as AcquRoot 
  class TMCGenerator (J.R.M.Annand -- http://nuclear.gla.ac.uk/~acqusys/doc).

  -- C.Tarbert
*/

#include <vector>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TMath.h"

#include "AMPTOOLS_MCGEN/NBodyPhaseSpaceFactory.h"

const double NBodyPhaseSpaceFactory::kPi = 3.14159;

NBodyPhaseSpaceFactory::NBodyPhaseSpaceFactory( double parentMass, const vector<double>& childMass ) :
  m_parentMass( parentMass ),
  m_childMass( childMass )
{
  m_Nd = (int)childMass.size();
}

vector<TLorentzVector>
NBodyPhaseSpaceFactory::generateDecay(bool uniformWeights) {
	

  vector<TLorentzVector> child( m_Nd );

  double Tcm = m_parentMass;
  int n, m;
  for( n=0; n<m_Nd; n++ ){ Tcm -= m_childMass[n]; }
  assert( Tcm > 0. );

  double emmax = Tcm + m_childMass[0];
  double emmin = 0;
  double wt = 1;
  for (n=1; n<m_Nd; n++) {
    emmin += m_childMass[n-1];
    emmax += m_childMass[n];
    wt *= pdk(emmax, emmin, m_childMass[n]);
  }
  double WtMax = 1/wt;                               // max weight for gating events

  double rnd[m_Nd];
  double invMas[m_Nd];
  double pd[m_Nd];
  int irnd[m_Nd];
  rnd[0] = 0.0; rnd[m_Nd-1] = 1.0;
  do{
    switch( m_Nd ){
    default:
      for( n=1; n<m_Nd-1; n++){ rnd[n] = random( 0., 1. ); }
      for( n=0; n<m_Nd; n++ ){ irnd[n] = n; }
      sort( irnd, irnd+m_Nd, CompareAsc<const double*>(rnd) ); // sort random numbers ascending
      rnd[m_Nd-1] = 1.0;
      break;
    case 3:                                     // 3 decay particles
      rnd[1] = random(0.0,1.0);
      irnd[0] = 0; irnd[1] = 1; irnd[2] = 2;
      break;
    case 2:
      irnd[0] = 0; irnd[1] = 1;
      break;
    }
    wt = 0.0;
    for (n=0; n<m_Nd; n++) {
      wt += m_childMass[n];
      invMas[n] = rnd[irnd[n]]*Tcm + wt;
    }
    wt = WtMax;
    for (n=0; n<m_Nd-1; n++) {
      pd[n] = pdk(invMas[n+1],invMas[n],m_childMass[n+1]);
      wt *= pd[n];
    }
  }while( uniformWeights && (random(0.0,WtMax) > wt) );

  if(uniformWeights) m_lastWt = 1.0;
  else               m_lastWt = wt;

  //
  // Specification of 4-momenta (Raubold-Lynch method)
  //
  child[0].SetPxPyPzE(0, pd[0], 0 , sqrt(pd[0]*pd[0]+m_childMass[0]*m_childMass[0]) );
  for(n=1;;){
    child[n].SetPxPyPzE(0, -pd[n-1], 0 ,
		    sqrt(pd[n-1]*pd[n-1]+m_childMass[n]*m_childMass[n]) );

    double cosZ = random(-1.,1.);
    double angY = random(0.0, 2.*kPi);
    for (m=0; m<=n; m++) {
      child[m].RotateZ( acos(cosZ) );
      child[m].RotateY( angY );
    }
    if( n == m_Nd-1 ) break;
    double beta = pd[n] / sqrt(pd[n]*pd[n] + invMas[n]*invMas[n]);
    for (m=0; m<=n; m++) child[m].Boost(0,beta,0);
    n++;
  }

  return child;
}

double
NBodyPhaseSpaceFactory::pdk( double a, double b, double c ) const {
	
  double x = (a-b-c)*(a+b+c)*(a-b+c)*(a+b-c);
  x = sqrt(x)/(2*a);
  return x;

}

double
NBodyPhaseSpaceFactory::random( double low, double hi ) const {
	
  return( ( hi - low ) * drand48() + low );
}
