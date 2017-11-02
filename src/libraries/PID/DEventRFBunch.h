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

		void Set_Content(DetectorSystem_t locTimeSource, double locTime, double locTimeVariance, unsigned int locNumParticleVotes = 0);
		void Reset(void);

		DetectorSystem_t dTimeSource = SYS_NULL; //e.g. SYS_RF, SYS_START, SYS_NULL (no valid source or not enough tracks/showers to pick it)

		double dTime; //The RF time propagated to the center of the target.  This time is defined at the center of the target. 
		double dTimeVariance;
		unsigned int dNumParticleVotes; //e.g. will trust time much more if 2+ rather than "1" or "0"

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "Source System", "%s", SystemName(dTimeSource));
			AddString(items, "t", "%3.5f", dTime);
			AddString(items, "var_t", "%3.2f", dTimeVariance);
			AddString(items, "#tracks", "%i", dNumParticleVotes);
		}
};

inline void DEventRFBunch::Reset(void)
{
	dTimeSource = SYS_NULL;
	dTime = 0.0;
	dTimeVariance = 0.0;
	dNumParticleVotes = 0;
}

inline void DEventRFBunch::Set_Content(DetectorSystem_t locTimeSource, double locTime, double locTimeVariance, unsigned int locNumParticleVotes)
{
	dTimeSource = locTimeSource;
	dTime = locTime;
	dTimeVariance = locTimeVariance;
	dNumParticleVotes = locNumParticleVotes;
}

#endif // _DEventRFBunch_

