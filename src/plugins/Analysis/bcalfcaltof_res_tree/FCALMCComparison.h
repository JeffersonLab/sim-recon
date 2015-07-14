//
//    File: FCALMCComparison.h
// Created: Thu Oct 29 09:49:51 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _FCALMCComparison_
#define _FCALMCComparison_

#include <TObject.h>

class FCALMCComparison : public TObject{

	public:

		FCALMCComparison(){};
		~FCALMCComparison(){};

		// Data members
		float dTrueX;
		float dTrueY;
		float dTrueZ;
		float dTrueE;
		float dTrueT;

		float dDeltaX;
		float dDeltaY;
		float dDeltaZ;
		float dDeltaE;
		float dDeltaT;

		float dShowerUncertaintyX;
		float dShowerUncertaintyY;
		float dShowerUncertaintyZ;
		float dShowerUncertaintyT;
		float dShowerUncertaintyE;

		float dPathLengthCorrection;

	private:
		ClassDef(FCALMCComparison, 1);

};

#endif // _FCALMCComparison_
