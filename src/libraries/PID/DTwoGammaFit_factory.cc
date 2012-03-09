// $Id: DTwoGammaFit_factory.cc 2496 2007-03-12 00:35:46Z kornicer $
//
//    File: DTwoGammaFit_factory.cc
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#include <math.h>

#include <JANA/JEvent.h>
using namespace jana;

#include "DTwoGammaFit.h"
#include "DKinFit.h"
#include "DVertex.h"
#include "DNeutralParticleHypothesis.h"
#include "DTwoGammaFit_factory.h"


bool SortByProb(const DTwoGammaFit *a,const DTwoGammaFit *b){
  return (a->getProb() > b->getProb());
}


//----------------
// Constructor
//----------------
DTwoGammaFit_factory::DTwoGammaFit_factory(double aMass)
{
  // Set defaults
  if (aMass > 0) { 
    fMass = aMass; 
  }
  else {
    cout << "Mass undefined, set to pi0" << endl; 
    fMass = 0.135; 
  }

}


//------------------
// evnt
// TwoGammaGit factory: loop over all pair combinations and fit them to the supplied mass using
//			kinemati fitter. Tags, like PI0, ETA (and ETAP eventually) can be used in derived 
//			factories.
//------------------
jerror_t DTwoGammaFit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{ 
  //so the way this works is that for a given DBCALShower/DFCALShower, a
  //DNeutralParticleHypothesis is created for each DVertex (DVertex's are
  //determined by grouping together DChargedTrackHypothesis). We
  //obviously don't want to consider combinations of hypotheses from
  //different vertices.

  vector<const DVertex*> vertices;
  vector<const DNeutralParticleHypothesis*> neutrals;
  eventLoop->Get(vertices);
  eventLoop->Get(neutrals);
   
  JObject::oid_t nPairs=0;
  //Loop over all fit candidates (all combinations of two
  //photon-hypotheses from a common vertex)
  for (unsigned int h=0; h<vertices.size(); h++) {
    for (unsigned int i=0; i<neutrals.size(); i++) {

      //only consider photon hypotheses
      if ( neutrals[i]->dPID != Gamma ) continue;

      //each hypothesis has a vertex associated with it
      vector<const DVertex*> hypothesisVertex;
      neutrals[i]->Get(hypothesisVertex);
      //each hypothesis should have exactly one vertex asscoiated with it
      //if not, something's wrong
      if (hypothesisVertex.size() != 1) continue;

      //check if this hypothesis comes from the vertex currently under
      //consideration
      if (hypothesisVertex[0] != vertices[h]) continue;

      for (unsigned int j = i+1; j < neutrals.size() ; j++) {
	//same checks as above
	if ( neutrals[j]->dPID != Gamma ) continue;
	neutrals[j]->Get(hypothesisVertex);
	if (hypothesisVertex.size() != 1) continue;
	if (hypothesisVertex[0] != vertices[h]) continue;

	nPairs++;

	//set up the fit
	DKinFit *kfit = new DKinFit(); 
	kfit->SetVerbose(0);

	vector<DKinematicData> kindata;
	kindata.clear();
	kindata.push_back( *(neutrals[i]) );
	kindata.push_back( *(neutrals[j]) );

	kfit->SetFinal(kindata);
	//fit!
	kfit->FitTwoGammas(fMass,1.); // second parameter is scale for photon error matrix
 
	vector<DKinematicData> kinout = kfit->GetFinal_out();
	if (kinout.size() != 2) continue;

	DTwoGammaFit* fit2g = new DTwoGammaFit( nPairs );

	// set two gamma kinematics
	DLorentzVector P4U = kindata[0].lorentzMomentum() + kindata[1].lorentzMomentum();
	DLorentzVector P4 = kinout[0].lorentzMomentum() + kinout[1].lorentzMomentum();
	DVector3 vertex = vertices[h]->dSpacetimeVertex.Vect();

	DMatrixDSym errorMatrix(7);

	for (int k=0; k<3; k++) {
	  for (int l=0; l<3; l++) {
	    errorMatrix[k][l]=kinout[0].errorMatrix()[k][l] + kinout[1].errorMatrix()[k][l];
	  }
	}

	fit2g->setUMass( P4U.M() );
	fit2g->setMass( P4.M() );
	fit2g->setMomentum( P4.Vect() );
	fit2g->setErrorMatrix(errorMatrix);
	fit2g->setPosition( vertex );
	fit2g->setProb( kfit->Prob());
	fit2g->setChi2( kfit->Chi2());
	fit2g->setNdf( kfit->Ndf());
	for (int k=0; k < 6; k++) {
	  fit2g->setPulls( kfit->GetPull(k),  k);
	}
	fit2g->setChildID( kindata[0].id, 0);
	fit2g->setChildID( kindata[1].id, 1);
	fit2g->setChildMom( kindata[0].lorentzMomentum() , 0);
	fit2g->setChildMom( kindata[1].lorentzMomentum() , 1);

	fit2g->setChildFit( kinout[0] , 0);
	fit2g->setChildFit( kinout[1] , 1);

	//we add the DVertex as an associated object to ease comparison of pi0
	//vertices with charged track vertices (you need only compare DVertex
	//pointers, rather than doing a fuzzy comparison of vertex coordinates)
	fit2g->AddAssociatedObject(vertices[h]);
 
	//must check if nan to prevent problems (crashes) in sort()
	//still need to look into avoiding nan's in first place
	if ( !isnan(fit2g->getProb()) ) _data.push_back( fit2g );
 
	delete kfit;

      } 
    }
  }

  // sort by probability so most probable fits are listed first
  sort(_data.begin(), _data.end(), SortByProb);

  return NOERROR;
}
