/*
 *  GammaPToXYP.cc
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "TLorentzVector.h"

#include "AMPTOOLS_MCGEN/GammaPToXYP.h"
#include "AMPTOOLS_MCGEN/TwoBodyDecayFactory.h"

#include "IUAmpTools/Kinematics.h"

// FORTRAN routines
extern "C"{
void cobrems_(float* Emax, float* Epeak, float* Dist);
float dntdx_(float* x);
};

// Wrapper function for total
double dNtdx(double x)
{
        float xx = x;
        return (double)dntdx_(&xx);
}

GammaPToXYP::GammaPToXYP( float lowMassXY, float highMassXY, 
                          float massX, float massY, float beamE,
                          ProductionMechanism::Type type ) : 
m_prodMech( ProductionMechanism::kProton, type, 6.0 ), // last arg is t dependence
m_beam( 0, 0, beamE, beamE),
m_target( 0, 0, 0, 0.938 ),
m_childMass( 0 ) {

  m_childMass.push_back( massX );
  m_childMass.push_back( massY );
  
  m_prodMech.setMassRange( lowMassXY, highMassXY );
 
  // Initialize coherent brem table
  float Emax =  5.5;
  float Epeak = 3.0;
  float Elow = 2.5;
  float Ehigh = 3.0;
  float Dist = 76.0;
  cobrems_(&Emax, &Epeak, &Dist);

  // Create histogram
  cobrem_vs_E = new TH1D("cobrem_vs_E", "Coherent Bremstrahlung vs. E_{#gamma}", 550, 0.0, 5.5);
  
  // Fill histogram
  for(int i=1; i<=cobrem_vs_E->GetNbinsX(); i++){
	  double x = cobrem_vs_E->GetBinCenter(i)/Emax;
	  double y = dNtdx(x);
	  if(cobrem_vs_E->GetBinCenter(i) > Elow && cobrem_vs_E->GetBinCenter(i) < Ehigh) 
		  cobrem_vs_E->SetBinContent(i, y);
  }

}

Kinematics* 
GammaPToXYP::generate(){

  double beamE = cobrem_vs_E->GetRandom();
  m_beam.SetPxPyPzE(0,0,beamE,beamE);

  TLorentzVector resonance = m_prodMech.produceResonance( m_beam );
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
GammaPToXYP::addResonance( float mass, float width, float bf ){
  
  m_prodMech.addResonance( mass, width, bf );
}

