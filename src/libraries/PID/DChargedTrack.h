// $Id$
//
//    File: DChargedTrack.h
// Created: Mon Dec  7 14:29:24 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrack_
#define _DChargedTrack_

#include <vector>

#include "TMath.h"

#include <JANA/JObject.h>

#include <particleType.h>
#include <PID/DChargedTrackHypothesis.h>

using namespace std;

class DChargedTrack : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DChargedTrack);

		oid_t candidateid; // unique id for this track //same as for DTrackCandidate
		vector<const DChargedTrackHypothesis*> dChargedTrackHypotheses;

		int Get_Charge(void) const;
		bool Contains_Charge(int locCharge) const;

		const DChargedTrackHypothesis* Get_Hypothesis(Particle_t locPID) const;
		const DChargedTrackHypothesis* Get_BestFOM(void) const;
		const DChargedTrackHypothesis* Get_BestTrackingFOM(void) const;

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "Nhypotheses", "%d", dChargedTrackHypotheses.size());
		}
};

inline bool DChargedTrack::Contains_Charge(int locCharge) const
{
	for(unsigned int loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i)
	{
		if(ParticleCharge(dChargedTrackHypotheses[loc_i]->PID()) == locCharge)
			return true;
	}
	return false;
}

inline int DChargedTrack::Get_Charge(void) const
{
	const DChargedTrackHypothesis* locChargedTrackHypothesis = Get_BestFOM();
	return ((locChargedTrackHypothesis == NULL) ? 0 : ParticleCharge(locChargedTrackHypothesis->PID()));
}

inline const DChargedTrackHypothesis* DChargedTrack::Get_Hypothesis(Particle_t locPID) const
{
	for(unsigned int loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i)
	{
		if(dChargedTrackHypotheses[loc_i]->PID() == locPID)
			return dChargedTrackHypotheses[loc_i];
	}
	return NULL;
}

inline const DChargedTrackHypothesis* DChargedTrack::Get_BestFOM(void) const
{
	if(dChargedTrackHypotheses.empty())
		return NULL;
	double locBestFOM = -2.0;
	const DChargedTrackHypothesis* locBestChargedTrackHypothesis = dChargedTrackHypotheses[0];
	for(size_t loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i)
	{
		if(dChargedTrackHypotheses[loc_i]->dFOM > locBestFOM)
		{
			locBestChargedTrackHypothesis = dChargedTrackHypotheses[loc_i];
			locBestFOM = locBestChargedTrackHypothesis->dFOM;
		}
	}
	return locBestChargedTrackHypothesis;
}

inline const DChargedTrackHypothesis* DChargedTrack::Get_BestTrackingFOM(void) const
{
	if(dChargedTrackHypotheses.empty())
		return NULL;
	double locBestFOM = -2.0;
	const DChargedTrackHypothesis* locBestChargedTrackHypothesis = dChargedTrackHypotheses[0];
	for(size_t loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i)
	{
		unsigned int locNDF = dChargedTrackHypotheses[loc_i]->dNDF_Track;
		double locFOM = (locNDF > 0) ? TMath::Prob(dChargedTrackHypotheses[loc_i]->dChiSq_Track, locNDF) : numeric_limits<double>::quiet_NaN();
		if(locFOM > locBestFOM)
		{
			locBestChargedTrackHypothesis = dChargedTrackHypotheses[loc_i];
			locBestFOM = locFOM;
		}
	}
	return locBestChargedTrackHypothesis;
}

#endif // _DChargedTrack_

