#ifndef _DKinFitResults_factory_
#define _DKinFitResults_factory_

#include <deque>
#include <vector>
#include <map>
#include <set>

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"

#include "KINFITTER/DKinFitter.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisResults.h"
#include "ANALYSIS/DKinFitResults.h"
#include "ANALYSIS/DKinFitUtils_GlueX.h"

using namespace std;
using namespace jana;

class DKinFitResults_factory : public jana::JFactory<DKinFitResults>
{
	public:
		DKinFitResults_factory(){};
		~DKinFitResults_factory(){};

	private:

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		size_t Get_NumKinFitReactions(JEventLoop* locEventLoop);
		DKinFitResults* Build_KinFitResults(const DParticleCombo* locParticleCombo, const DKinFitChain* locKinFitChain);

		DKinFitter* dKinFitter;
		DKinFitUtils_GlueX* dKinFitUtils;

		unsigned int dKinFitDebugLevel;

		map<set<DKinFitConstraint*>, DKinFitResults*> dConstraintResultsMap; //used for determining if kinfit results will be identical
};

#endif // _DKinFitResults_factory_

