// $Id$
//
//    File: DTrackHit.h
// Created: Tue Aug 23 05:00:03 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackHit_
#define _DTrackHit_

#include <cmath>

#include <DMatrix.h>

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

#include "GlueX.h"

class DTrackHit:public JObject{
	public:
		JOBJECT_PUBLIC(DTrackHit);
		
		void InitCovarianceMatrix(void);

		float x,y,z,r,phi;
		DetectorSystem_t system;
		DMatrix cov; // covariance matrix rotated into lab x,y,z

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "x(cm)", "%3.1f", x);
			AddString(items, "y(cm)", "%3.1f", y);
			AddString(items, "z(cm)", "%3.1f", z);
			AddString(items, "r(cm)", "%3.1f", r);
			AddString(items, "phi(deg)", "%3.1f", phi*180.0/M_PI);
			AddString(items, "system", "%s", SystemName(system));
		}
};

#endif // _DTrackHit_

