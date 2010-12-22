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
#include "IUAmpTools/AmpVecs.h"
#include "CLHEP/Vector/LorentzVector.h"

#ifdef	GPU_ACCELERATION
#include "GPUManager/GPUManager.h"
#endif	//GPU_ACCELERATION

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
GammaPToXYP::generateOne(){
  
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

AmpVecs*
GammaPToXYP::generateMany( int nEvents ){
  
  AmpVecs* a = new AmpVecs;
  
  a->m_iNTrueEvents = nEvents;
  a->m_iNParticles = 4;
  
  Kinematics* pKinEvent;
  
#ifndef GPU_ACCELERATION
	
	a->m_iNEvents = a->m_iNTrueEvents;	
  
	a->m_pdData = new GDouble[4*a->m_iNParticles*a->m_iNEvents];
	a->m_pdWeights = new GDouble[a->m_iNEvents];
	
	int iEvent, iParticle, iDataIndex=0;	
	for( iEvent=0; iEvent < a->m_iNEvents; iEvent++ )
	{	
		pKinEvent = generateOne();
    
		for( iParticle = 0; iParticle < a->m_iNParticles; iParticle++ )
		{
			a->m_pdData[iDataIndex++] = pKinEvent->particle(iParticle).e();
			a->m_pdData[iDataIndex++] = pKinEvent->particle(iParticle).px();
			a->m_pdData[iDataIndex++] = pKinEvent->particle(iParticle).py();
			a->m_pdData[iDataIndex++] = pKinEvent->particle(iParticle).pz();
		}
		
		a->m_pdWeights[iEvent] = pKinEvent->weight();
		
		delete pKinEvent;
	}
	
#else
  
	// GPU needs all the 
	// Padding the arrays to be suitable for GPU	I/O
	a->m_iNEvents = GPUManager::calcNEventsGPU(a->m_iNTrueEvents);		
  
	a->m_pdData = new GDouble[4*a->m_iNParticles*a->m_iNEvents];
	a->m_pdWeights = new GDouble[a->m_iNEvents];
	
	int iEvent, iParticle;	
	for( iEvent = 0; iEvent < a->m_iNEvents; iEvent++ )
	{
		// Skip the first read and after the last data event
		// Pad the rest of the array by copying the last event
		if( iEvent < a->m_iNTrueEvents )
			pKinEvent = generateOne();
    
		for( iParticle = 0; iParticle < a->m_iNParticles; iParticle++ )
		{
			a->m_pdData[(4*iParticle+0)*a->m_iNEvents+iEvent] = 
        pKinEvent->particle(iParticle).e();
			a->m_pdData[(4*iParticle+1)*a->m_iNEvents+iEvent] = 
        pKinEvent->particle(iParticle).px();
			a->m_pdData[(4*iParticle+2)*a->m_iNEvents+iEvent] = 
        pKinEvent->particle(iParticle).py();
			a->m_pdData[(4*iParticle+3)*a->m_iNEvents+iEvent] = 
        pKinEvent->particle(iParticle).pz();
		}
		
		a->m_pdWeights[iEvent] = pKinEvent->weight();
		
		//Must free the pointer if we're not at the last event
		if( iEvent < ( a->m_iNTrueEvents - 1 ) ) 
      delete pKinEvent;
	}
  
  // clean up the pointer for the last event if we've
  // used it to pad the array
  if( a->m_iNEvents > a->m_iNTrueEvents ) delete pKinEvent;
  
#endif // GPU_ACCELERATION	
  
  return a;
}

void
GammaPToXYP::addResonance( float mass, float width, float bf ){
  
  m_prodMech.addResonance( mass, width, bf );
}

