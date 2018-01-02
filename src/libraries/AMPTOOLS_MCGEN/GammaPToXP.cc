/*
 *  GammaPToXP.cc
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "TLorentzVector.h"
#include "TRandom3.h"

#include "AMPTOOLS_MCGEN/GammaPToXP.h"
#include "AMPTOOLS_MCGEN/BeamProperties.h"

GammaPToXP::GammaPToXP( float massX, TString beamConfigFile) : 
m_target( 0, 0, 0, 0.938 ),
m_childMass( 0 ) {

  m_childMass.push_back( massX );

  // get beam properties from configuration file
  BeamProperties beamProp(beamConfigFile);
  cobrem_vs_E = (TH1D*)beamProp.GetFlux();
  cobrem_vs_E->GetName();
}

Kinematics* 
GammaPToXP::generate(){

  double beamE = cobrem_vs_E->GetRandom();
  m_beam.SetPxPyPzE(0,0,beamE,beamE);
  TLorentzVector cm = m_beam + m_target;

  Double_t masses[2] = {0.938,m_childMass[0]};
  TGenPhaseSpace phsp;
  phsp.SetDecay(cm,2,masses);

  double phsp_wt_max = phsp.GetWtMax();
  double genWeight;
  do {
     genWeight = phsp.Generate();
  }
  while( random(0., phsp_wt_max) >= genWeight || genWeight != genWeight);
  
  TLorentzVector *recoil = phsp.GetDecay(0);
  TLorentzVector *pX = phsp.GetDecay(1);

  vector< TLorentzVector > allPart;
  allPart.push_back( m_beam );
  allPart.push_back( *recoil );
  allPart.push_back( *pX );
 
  return new Kinematics( allPart, genWeight );
}

double
GammaPToXP::random( double low, double hi ) const {

        return( ( hi - low ) * gRandom->Uniform() + low );
}

