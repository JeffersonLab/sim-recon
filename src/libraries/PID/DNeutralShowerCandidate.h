// $Id$
//
//    File: DNeutralShowerCandidate.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralShowerCandidate_
#define _DNeutralShowerCandidate_

#include <vector>
#include <JANA/JObject.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <GlueX.h>
#include <DLorentzVector.h>

using namespace std;

class DNeutralShowerCandidate : public jana::JObject{
	public:
		JOBJECT_PUBLIC(DNeutralShowerCandidate);

		DLorentzVector dSpacetimeVertex;
		DLorentzVector dSpacetimeVertexUncertainties;
		float dEnergy;
		float dEnergyUncertainty;
		DetectorSystem_t dDetectorSystem;

		DNeutralShowerCandidate(const DBCALShower *locBCALShower){
			dDetectorSystem = SYS_BCAL;
			dEnergy = locBCALShower->E;
			dEnergyUncertainty = (dEnergy >= 0.0) ? pow( 0.0445*sqrt( dEnergy ) + 0.009*dEnergy, 2.0) : 1e-6; //from old DPhoton_factory::makeBCalPhoton() function
			dSpacetimeVertex.SetXYZT(locBCALShower->x, locBCALShower->y, locBCALShower->z, locBCALShower->t);
			dSpacetimeVertexUncertainties.SetXYZT(locBCALShower->xErr, locBCALShower->yErr, locBCALShower->zErr, locBCALShower->tErr);
			this->AddAssociatedObject(locBCALShower);
		}

		DNeutralShowerCandidate(const DFCALShower *locFCALShower){
			dDetectorSystem = SYS_FCAL;
			dEnergy = locFCALShower->getEnergy();
			dEnergyUncertainty = (dEnergy >= 0.0) ? pow( 0.042*sqrt(dEnergy) + 0.0001, 2.0 ) : 1e-6; //from old DPhoton_factory::makeFCalPhoton() function
			dSpacetimeVertex.SetVect(locFCALShower->getPosition());
			dSpacetimeVertex.SetT(locFCALShower->getTime());
			dSpacetimeVertexUncertainties.SetVect(locFCALShower->getPositionError());
			dSpacetimeVertexUncertainties.SetT(0.0); //not stored in DFCALShower
			this->AddAssociatedObject(locFCALShower);
		}
};

#endif // _DNeutralShowerCandidate_

