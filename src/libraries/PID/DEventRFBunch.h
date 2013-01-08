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

#include "JANA/JObject.h"

using namespace std;

class DEventRFBunch : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DEventRFBunch);

		double dTime; //The RF time propagated to the center of the target.  This time is defined at the center of the target. 
		double dTimeVariance;
		bool dMatchedToTracksFlag; //true if confident in value from matching to tracks, otherwise false (just keeps raw value)

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "t", "%3.5f", dTime);
			AddString(items, "var_t", "%3.2f", dTimeVariance);
		}
};

#endif // _DEventRFBunch_

