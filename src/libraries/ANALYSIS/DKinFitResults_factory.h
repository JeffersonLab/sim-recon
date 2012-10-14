#ifndef _DKinFitResults_factory_
#define _DKinFitResults_factory_

#include <deque>
#include <vector>
#include <map>

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"

#include "HDGEOMETRY/DMagneticFieldMap.h"

#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "PID/DKinematicData.h"
#include "PID/DBeamPhoton.h"

#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DKinFitResults.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DKinFitter_GlueX.h"
#include "ANALYSIS/DAnalysisResults.h"

using namespace std;
using namespace jana;

class DKinFitResults_factory : public jana::JFactory<DKinFitResults>
{
	public:
		DKinFitResults_factory(){};
		~DKinFitResults_factory(){};

		void Reset_NewEvent(void);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void Setup_KinFit(const DParticleCombo* locParticleCombo, deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles);
		void Setup_P4Constraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_P4, deque<const DKinFitParticle*>& locFinalKinFitParticles_P4, deque<size_t>& locIncludedStepIndices);
		void Setup_VertexConstraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_Vertex, deque<const DKinFitParticle*>& locFinalKinFitParticles_Vertex, deque<size_t>& locIncludedStepIndices);

		TVector3 Calc_VertexGuess(const deque<const DKinFitParticle*>& locInitialKinFitParticles, const deque<const DKinFitParticle*>& locFinalKinFitParticles);
		double Calc_TimeGuess(const deque<const DKinFitParticle*>& locFinalKinFitParticles, DVector3 locVertexGuess, bool locUseRFTimeFlag, double locRFTime);
		void Remove_BadVertexConstraints(deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_Vertices, deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_Vertices) const;

		void Build_KinFitResults(const DParticleCombo* locParticleCombo, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_Input, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_Input);

		const DAnalysisUtilities* dAnalysisUtilities;
		DKinFitter_GlueX dKinFitter;

		unsigned int dDebugLevel;
		unsigned int dKinFitDebugLevel;
		bool dLinkVerticesFlag;
};

#endif // _DKinFitResults_factory_

