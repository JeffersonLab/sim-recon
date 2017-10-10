#ifndef _DChargedTrack_factory_Combo_
#define _DChargedTrack_factory_Combo_

#include <string>
#include <vector>
#include <unordered_map>

#include <JANA/JFactory.h>
#include <TRACKING/DTrackTimeBased.h>
#include "PID/DEventRFBunch.h"
#include "PID/DDetectorMatches.h"
#include <PID/DChargedTrack.h>
#include <PID/DChargedTrackHypothesis.h>
#include "PID/DChargedTrackHypothesis_factory.h"

using namespace std;
using namespace jana;

class DChargedTrack_factory_Combo : public jana::JFactory<DChargedTrack>
{
	public:
		DChargedTrack_factory_Combo(){};
		~DChargedTrack_factory_Combo(){};
		const char* Tag(void){return "Combo";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.

		string dTrackSelectionTag;
		DChargedTrackHypothesis_factory* dChargedTrackHypothesisFactory;
		vector<DChargedTrackHypothesis*> dCreatedHypotheses;
};

#endif // _DChargedTrack_factory_Combo_

