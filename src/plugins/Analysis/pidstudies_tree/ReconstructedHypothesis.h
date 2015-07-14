//
//    File: ReconstructedHypothesis.h
// Created: Thu Oct 29 09:49:51 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _ReconstructedHypothesis_
#define _ReconstructedHypothesis_

#include <TObject.h>
#include <particleType.h>
#include <DLorentzVector.h>

class ReconstructedHypothesis : public TObject{

	public:

		ReconstructedHypothesis(){};
		~ReconstructedHypothesis(){};

		// Data members
		Particle_t dPID;
		DLorentzVector dFourMomentum;
		DLorentzVector dSpacetimeVertex; // vertex position in cm + vertex time in ns

		double dChiSq_Overall;
		int dNDF_Overall;

		double dChiSq_Tracking;
		int dNDF_Tracking;

		double dChiSq_DCdEdx;
		int dNDF_DCdEdx;

		double dChiSq_Timing;
		int dNDF_Timing;

		double dChiSq_TOFdEdx;
		int dNDF_TOFdEdx;

		double dChiSq_BCALdEdx;
		int dNDF_BCALdEdx;

	private:
		ClassDef(ReconstructedHypothesis, 1);

};

#endif // _ReconstructedHypothesis_
