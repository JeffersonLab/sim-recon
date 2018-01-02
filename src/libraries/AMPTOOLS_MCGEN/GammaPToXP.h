#if !defined(GAMMAPTOXP)
#define GAMMAPTOXP

/*
 *  GammaPToXP.h
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH1.h"

#include "IUAmpTools/Kinematics.h"

class Kinematics;

class GammaPToXP {
  
public:
  
  GammaPToXP( float massX, TString beamConfigFile);
  
  Kinematics* generate();
  
private:
  
  TLorentzVector m_beam;
  TLorentzVector m_target;
  
  vector< double > m_childMass;
  
  TH1D *cobrem_vs_E;
  
  double random( double low, double hi ) const;

};

#endif
