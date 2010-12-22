#if !defined(GAMMAPTOXYP)
#define GAMMAPTOXYP

/*
 *  GammaPToXYP.h
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "CLHEP/Vector/LorentzVector.h"

class Kinematics;
class AmpVecs;

class GammaPToXYP {
  
public:
  
  GammaPToXYP( float lowMassXY, float highMassXY, float massX, float massY,
               ProductionMechanism::Type type );
  
  Kinematics* generateOne();
  AmpVecs* generateMany( int nEvents );
  
  void addResonance( float mass, float width, float bf );
  
private:

  ProductionMechanism m_prodMech;
  
  HepLorentzVector m_beam;
  HepLorentzVector m_target;
  
  vector< double > m_childMass;
  
};

#endif
