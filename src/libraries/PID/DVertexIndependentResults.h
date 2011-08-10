// $Id$
//
//    File: DVertexIndependentResults.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DVertexIndependentResults_
#define _DVertexIndependentResults_

#include <vector>
#include <JANA/JObject.h>
#include <PID/DNeutralShowerCandidate.h>
#include <PID/DChargedTrack.h>

using namespace std;

class DVertexIndependentResults : public jana::JObject {
	public:
		JOBJECT_PUBLIC(DVertexIndependentResults);

		vector <const DChargedTrack*> dChargedTracks;
		vector <const DNeutralShowerCandidate*> dNeutralShowerCandidates;
};

#endif // _DVertexIndependentResults_

