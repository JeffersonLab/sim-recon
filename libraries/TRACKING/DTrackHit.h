// $Id$
//
//    File: DTrackHit.h
// Created: Tue Aug 23 05:00:03 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackHit_
#define _DTrackHit_

#include <DMatrix.h>

#include "JANA/JObject.h"
#include "JANA/JFactory.h"
#include "GlueX.h"

class DTrackHit:public JObject{
	public:
		HDCLASSDEF(DTrackHit);
		
		void InitCovarianceMatrix(void);

		float x,y,z,r,phi;
		DetectorSystem_t system;
		DMatrix cov; // covariance matrix rotated into lab x,y,z
};

#endif // _DTrackHit_

