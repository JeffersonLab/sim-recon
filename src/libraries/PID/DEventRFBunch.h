// $Id$
//
//    File: DEventRFBunch.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DEventRFBunch_
#define _DEventRFBunch_

#include <vector>
#include <utility>
#include <string>

#include "GlueX.h"
#include "JANA/JObject.h"

using namespace std;

class DEventRFBunch : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DEventRFBunch);

		DetectorSystem_t dTimeSource; //e.g. SYS_RF, SYS_START

		double dTime; //The RF time propagated to the center of the target.  This time is defined at the center of the target. 
		double dTimeVariance;
		unsigned int dNumParticleVotes; //e.g. will trust time much more if 2+ rather than "1" or "0"

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "t", "%3.5f", dTime);
			AddString(items, "var_t", "%3.2f", dTimeVariance);
			AddString(items, "#tracks", "%i", dNumParticleVotes);
		}
};

#endif // _DEventRFBunch_

