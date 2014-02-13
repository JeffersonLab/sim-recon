/*
 *  GammaPToXYP.cc
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "AMPTOOLS_MCGEN/GammaPToXYP.h"
#include "AMPTOOLS_MCGEN/TwoBodyDecayFactory.h"

#include "IUAmpTools/Kinematics.h"
#include "CLHEP/Vector/LorentzVector.h"

GammaPToXYP::GammaPToXYP( float lowMassXY, float highMassXY, 
                          float massX, float massY,
                          ProductionMechanism::Type type ) : 
m_prodMech( ProductionMechanism::kProton, type, 4.0 ), // last arg is t dependence
m_beam( 0, 0, 9, 9 ),
m_target( 0, 0, 0, 0.938 ),
m_childMass( 0 ) {
  
  m_childMass.push_back( massX );
  m_childMass.push_back( massY );
  
  m_prodMech.setMassRange( lowMassXY, highMassXY );
  
}

Kinematics* 
GammaPToXYP::generate(){
  
  HepLorentzVector resonance = m_prodMech.produceResonance( m_beam );
  double genWeight = m_prodMech.getLastGeneratedWeight();
  
  vector< HepLorentzVector > allPart;
  allPart.push_back( m_beam );
  allPart.push_back( m_beam + m_target - resonance );
  
  TwoBodyDecayFactory decay( resonance.m(), m_childMass );
  
  vector<HepLorentzVector> fsPart = decay.generateDecay();
  
  for( vector<HepLorentzVector>::iterator aPart = fsPart.begin();
      aPart != fsPart.end(); ++aPart ){
    
    aPart->boost( resonance.boostVector() );
    allPart.push_back( *aPart );
  }
  
  return new Kinematics( allPart, genWeight );
}

void
GammaPToXYP::addResonance( float mass, float width, float bf ){
  
  m_prodMech.addResonance( mass, width, bf );
}

