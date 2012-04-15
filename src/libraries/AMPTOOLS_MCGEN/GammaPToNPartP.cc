/*
 *  GammaPToNPartP.h
 *   by Igor Senderovich
 *  structure based on GammaToXYZP
 *  written by Matthew Shepherd 
 */

#include "GammaPToNPartP.h"
#include "particleType.h"
#include "AMPTOOLS_MCGEN/DalitzDecayFactory.h"
#include "TGenPhaseSpace.h"
#include "NBodyPhaseSpaceFactory.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/AmpVecs.h"
#include "CLHEP/Vector/LorentzVector.h"

#ifdef	GPU_ACCELERATION
#include "GPUManager/GPUManager.h"
#endif	//GPU_ACCELERATION


GammaPToNPartP::GammaPToNPartP( float lowMass, float highMass, 
				vector<double> &ChildMass,
				ProductionMechanism::Type type, 
				float tcoef, float Ebeam) : 
  m_prodMech( ProductionMechanism::kProton, type, tcoef ), // last arg is t dependence
  m_beam( 0, 0, Ebeam, Ebeam ),
  m_target( 0, 0, 0, ParticleMass(Proton) ),
  m_ChildMass(ChildMass)
{
  assert(Ebeam>0);
  
  m_Npart = ChildMass.size();
  assert(m_Npart>0);

  m_prodMech.setMassRange( lowMass, highMass );


}

/**
 * The function generates a N particle final
 * state event consistent with N-body phase space.
 * (No intermediate resonances are used for important sampling.)
 */
Kinematics* 
GammaPToNPartP::generateOne(){

  HepLorentzVector resonance;
  do{
    resonance=m_prodMech.produceResonance( m_beam );
  }while(!(resonance.e() < m_beam.e()));


  //TLorentzVector tresonance(resonance.px(),resonance.py(),
  //		    resonance.pz(),resonance.e());
  double genWeight = m_prodMech.getLastGeneratedWeight();
  
  vector< HepLorentzVector > allPart;
  allPart.push_back( m_beam );
  allPart.push_back( m_beam + m_target - resonance );
  
  // X decay phase space 
  /*TGenPhaseSpace Xdecay;
  Xdecay.SetDecay(tresonance, m_Npart, m_ChildMass);
  genWeight *= Xdecay.Generate();
  */

  NBodyPhaseSpaceFactory psFactory(resonance.m(),m_ChildMass);
  vector< HepLorentzVector > children = psFactory.generateDecay(false);
  genWeight *= psFactory.getLastGeneratedWeight();

  Hep3Vector b3(resonance.boostVector());   // boost vector from parent
  for (unsigned int n=0; n<children.size(); ++n ){
    children[n].boost(b3);
    allPart.push_back(children[n]);
  }

  /*
  for(unsigned int i = 0 ; i<m_Npart ; ++i){
    TLorentzVector *tPart = Xdecay.GetDecay(i);
    HepLorentzVector Part(tPart->Px(),tPart->Py(),tPart->Pz(),tPart->Energy());
    allPart.push_back(Part);
    }*/
  
  return new Kinematics( allPart, genWeight );
}

AmpVecs*
GammaPToNPartP::generateMany( int nEvents ){
  
  AmpVecs* a = new AmpVecs;
  
  a->m_iNTrueEvents = nEvents;
  a->m_iNParticles = m_Npart+2;
  
  Kinematics* pKinEvent;
  
#ifndef GPU_ACCELERATION
  
  a->m_iNEvents = a->m_iNTrueEvents;	
  
  a->m_pdData = new GDouble[4*a->m_iNParticles*a->m_iNEvents];
  a->m_pdWeights = new GDouble[a->m_iNEvents];
  
  int iEvent, iParticle, iDataIndex=0;	
  for( iEvent=0; iEvent < a->m_iNEvents; iEvent++ ) {	
    pKinEvent = generateOne();
    
    for( iParticle = 0; iParticle < a->m_iNParticles; iParticle++ ) {
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
  for( iEvent = 0; iEvent < a->m_iNEvents; iEvent++ ){
    // Skip the first read and after the last data event
    // Pad the rest of the array by copying the last event
    if( iEvent < a->m_iNTrueEvents )
      pKinEvent = generateOne();
    
    for( iParticle = 0; iParticle < a->m_iNParticles; iParticle++ ){
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
GammaPToNPartP::addResonance( float mass, float width, float bf ){
  
  m_prodMech.addResonance( mass, width, bf );
}

