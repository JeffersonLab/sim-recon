#ifndef _DTrackTimeBased_factory_Combo_
#define _DTrackTimeBased_factory_Combo_

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "JANA/JFactory.h"
#include "particleType.h"

#include "TRACKING/DTrackTimeBased.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"

#include "ANALYSIS/DReaction.h"

using namespace jana;
using namespace std;

class DTrackTimeBased_factory_Combo:public jana::JFactory<DTrackTimeBased>
{
	public:
		DTrackTimeBased_factory_Combo(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DTrackTimeBased_factory_Combo(){};
		const char* Tag(void){return "Combo";}

		vector<Particle_t> Get_ParticleIDsToTry(Particle_t locPID) const
		{
			auto locIterator = dParticleIDsToTry.find(locPID);
			if(locIterator == dParticleIDsToTry.end())
				return vector<Particle_t>();
			return locIterator->second;
		}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.

		void Create_PIDsAsNeeded(const DChargedTrack* locChargedTrack, set<Particle_t>& locPIDs);
		const DChargedTrackHypothesis* Get_ChargedHypothesisToUse(const DChargedTrack* locChargedTrack, Particle_t locDesiredPID);
		DTrackTimeBased* Convert_ChargedTrack(const DChargedTrackHypothesis* locChargedTrackHypothesis, Particle_t locNewPID);

		map<Particle_t, vector<Particle_t> > dParticleIDsToTry;

		set<Particle_t> dPositivelyChargedPIDs;
		set<Particle_t> dNegativelyChargedPIDs;

		string dTrackSelectionTag;
};

#endif // _DTrackTimeBased_factory_Combo_
