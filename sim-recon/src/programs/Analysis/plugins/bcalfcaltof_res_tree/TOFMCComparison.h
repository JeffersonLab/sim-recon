//
//    File: TOFMCComparison.h
// Created: Thu Oct 29 09:49:51 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _TOFMCComparison_
#define _TOFMCComparison_

#include <TObject.h>

class TOFMCComparison : public TObject{

	public:

		TOFMCComparison(){};
		~TOFMCComparison(){};

		// Data members
		float dTrueX;
		float dTrueY;
		float dTrueZ;
		float dTruedE;
		float dTrueT;
		float dTrueBetaGamma;

		float dDeltaX;
		float dDeltaY;
		float dDeltaZ;
		float dDeltadE;
		float dDeltaT;

		float dPathLengthCorrection;
		bool dHorizontalPlaneFlag;
		bool dVerticalPlaneFlag;

	private:
		ClassDef(TOFMCComparison, 1);

};

#endif // _TOFMCComparison_
