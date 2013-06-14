/*
 *  GammaPToXYZP.cc
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 5/25/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "GammaPToXYZP.h"

#include "AMPTOOLS_MCGEN/DalitzDecayFactory.h"
#include "IUAmpTools/Kinematics.h"
#include "CLHEP/Vector/LorentzVector.h"

GammaPToXYZP::GammaPToXYZP( float lowMassXYZ, float highMassXYZ, 
                            float massX, float massY, float massZ,
                            ProductionMechanism::Type type ) : 
m_prodMech( ProductionMechanism::kProton, type, 4.0 ), // last arg is t dependence
m_beam( 0, 0, 9, 9 ),
m_target( 0, 0, 0, 0.938 ),
m_childMass( 0 ) {
  
  m_childMass.push_back( massX );
  m_childMass.push_back( massY );
  m_childMass.push_back( massZ );
  
  m_prodMech.setMassRange( lowMassXYZ, highMassXYZ );
  
}

Kinematics* 
GammaPToXYZP::generate(){
  
  HepLorentzVector resonance = m_prodMech.produceResonance( m_beam );
  double genWeight = m_prodMech.getLastGeneratedWeight();
  
  vector< HepLorentzVector > allPart;
  allPart.push_back( m_beam );
  allPart.push_back( m_beam + m_target - resonance );
  
  DalitzDecayFactory decay( resonance.m(), m_childMass );
  
  vector<HepLorentzVector> fsPart = decay.generateDecay();
  
  for( vector<HepLorentzVector>::iterator aPart = fsPart.begin();
      aPart != fsPart.end(); ++aPart ){
    
    aPart->boost( resonance.boostVector() );
    allPart.push_back( *aPart );
  }
  
  return new Kinematics( allPart, genWeight );
}

void
GammaPToXYZP::addResonance( float mass, float width, float bf ){
  
  m_prodMech.addResonance( mass, width, bf );
}

