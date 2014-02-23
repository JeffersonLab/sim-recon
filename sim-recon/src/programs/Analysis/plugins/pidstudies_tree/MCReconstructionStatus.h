//
//    File: MCReconstructionStatus.h
// Created: Thu Oct 29 09:49:51 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _MCReconstructionStatus_
#define _MCReconstructionStatus_

#include <TObject.h>
#include <vector>
#include <particleType.h>
#include <DLorentzVector.h>
#include <ReconstructedHypothesis.h>

class MCReconstructionStatus : public TObject{

	public:

		MCReconstructionStatus(){};
		~MCReconstructionStatus(){};

		// Data members
		Particle_t dThrownID;
		DLorentzVector dThrownFourMomentum;
		DLorentzVector dThrownSpacetimeVertex; // vertex position in cm + vertex time in ns

		vector<ReconstructedHypothesis*> dReconstructedHypothesisVector;

	private:
		ClassDef(MCReconstructionStatus, 1);

};

#endif // _MCReconstructionStatus_
