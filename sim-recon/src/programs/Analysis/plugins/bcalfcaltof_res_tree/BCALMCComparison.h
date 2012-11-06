//
//    File: BCALMCComparison.h
// Created: Thu Oct 29 09:49:51 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _BCALMCComparison_
#define _BCALMCComparison_

#include <TObject.h>

class BCALMCComparison : public TObject{

	public:

		BCALMCComparison(){};
		~BCALMCComparison(){};

		// Data members
		float dTrueR;
		float dTrueZ;
		float dTruePhi;
		float dTrueE;
		float dTrueT;

		float dDeltaR;
		float dDeltaZ;
		float dDeltaPhi;
		float dDeltaE;
		float dDeltaT;

		float dShowerUncertaintyX;
		float dShowerUncertaintyY;
		float dShowerUncertaintyZ;
		float dShowerUncertaintyT;
		float dShowerUncertaintyE;

		float dPathLengthCorrection;

	private:
		ClassDef(BCALMCComparison, 1);

};

#endif // _BCALMCComparison_
