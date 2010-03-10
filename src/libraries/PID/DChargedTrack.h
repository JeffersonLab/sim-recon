// $Id$
//
//    File: DChargedTrack.h
// Created: Mon Dec  7 14:29:24 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrack_
#define _DChargedTrack_

#include <vector>
#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <TRACKING/DTrackTimeBased.h>
using namespace std;

class DChargedTrack:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DChargedTrack);

		vector<const DTrackTimeBased *>hypotheses;
		                
	

};

#endif // _DChargedTrack_

