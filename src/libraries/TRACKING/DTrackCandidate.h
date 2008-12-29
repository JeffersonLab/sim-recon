// $Id$
//
//    File: DTrackCandidate.h
// Created: Sun Apr  3 12:38:16 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_
#define _DTrackCandidate_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

#include "PID/DKinematicData.h"

#define MAX_IHITS 256

class DTrackCandidate:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DTrackCandidate);

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "id", "0x%x", id);
		}
};

#endif // _DTrackCandidate_

