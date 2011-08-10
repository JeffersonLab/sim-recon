// $Id$
//
//    File: DNeutralTrack.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralTrack_
#define _DNeutralTrack_

#include <vector>
#include <JANA/JObject.h>
#include <PID/DNeutralTrackHypothesis.h>

using namespace std;

class DNeutralTrack:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DNeutralTrack);

		vector<const DNeutralTrackHypothesis *> dNeutralTrackHypotheses;

		void toStrings(vector<pair<string,string> > &items) const{
			AddString(items, "Nhypotheses", "%d", dNeutralTrackHypotheses.size());
		}

};

#endif // _DNeutralTrack_

