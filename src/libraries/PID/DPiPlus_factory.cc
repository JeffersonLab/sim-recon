// $Id$
//
//    File: DPiPlus_factory.cc
// Created: Sat Apr 14 12:13:05 EDT 2012
// Creator: davidl (on Darwin genmacbook.local 11.3.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <PID/DChargedTrack.h>

#include "DPiPlus_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPiPlus_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPiPlus_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPiPlus_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get list of all charged tracks
	vector<const DChargedTrack*> chargedtracks;
	loop->Get(chargedtracks);
	
	// Loop over charged tracks looking for ones where the most
	// probable hypothesis is a pi+
	for(unsigned int i=0; i<chargedtracks.size(); i++){
		const DChargedTrack *chargedtrack = chargedtracks[i];
		if(chargedtrack->dChargedTrackHypotheses.size()<1)continue;
		
		const DChargedTrackHypothesis* hypothesis = chargedtrack->dChargedTrackHypotheses[0];
		if(hypothesis->dPID == PiPlus){
			// Most probable hypothesis is a pi+. Make copy and store it in factory
			DPiPlus *pip = new DPiPlus(hypothesis);
			_data.push_back(pip);
		}
	}
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPiPlus_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPiPlus_factory::fini(void)
{
	return NOERROR;
}

