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

#include "TLorentzVector.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "TH1.h"

class Kinematics;
class AmpVecs;

class GammaPToXYP {
  
public:
  
  GammaPToXYP( float lowMassXY, float highMassXY, float massX, float massY,
               float beamE, ProductionMechanism::Type type );
  
  Kinematics* generate();
  
  void addResonance( float mass, float width, float bf );
  
private:

  ProductionMechanism m_prodMech;
  
  TLorentzVector m_beam;
  TLorentzVector m_target;
  
  vector< double > m_childMass;

  TH1D *cobrem_vs_E;
  
};

#endif
