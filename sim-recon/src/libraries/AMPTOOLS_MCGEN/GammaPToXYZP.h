#if !defined(GAMMAPTOXYZP)
#define GAMMAPTOXYZP

/*
 *  GammaPToXYZP.h
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 5/25/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */


#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "CLHEP/Vector/LorentzVector.h"

class Kinematics;
class AmpVecs;

class GammaPToXYZP {
  
public:
  
  GammaPToXYZP( float lowMassXYZ, float highMassXYZ, 
                float massX, float massY, float massZ,
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
