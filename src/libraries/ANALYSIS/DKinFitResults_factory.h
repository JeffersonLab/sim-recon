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

		bool Setup_KinFit(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locDecayingKinFitParticles);
		void Setup_P4Constraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_P4, deque<const DKinFitParticle*>& locFinalKinFitParticles_P4, deque<size_t>& locIncludedStepIndices);
		void Setup_VertexConstraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_Vertex, deque<const DKinFitParticle*>& locFinalKinFitParticles_Vertex, deque<size_t>& locIncludedStepIndices);

		bool Find_ConstrainableParticles(deque<pair<const DKinFitParticle*, size_t> >& locConstrainableParticles, const deque<const DKinFitParticle*>& locConstrainedParticles, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_P4s, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_P4s, const deque<size_t>& locConstraintsSetIndices);
		void Calc_P4Guess(pair<const DKinFitParticle*, size_t>& locConstrainableParticle, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_P4s, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_P4s, TVector3& locMomentum);

		bool Calc_VertexGuess(JEventLoop* locEventLoop, const deque<const DKinFitParticle*>& locInitialKinFitParticles, const deque<const DKinFitParticle*>& locFinalKinFitParticles, TVector3& locVertexGuess, deque<const DKinFitParticle*>& locFinalKinFitParticles_Vertex_Updated, int locVertexFindFlag, const map<const DKinFitParticle*, TVector3>& locP4GuessMap, map<const DKinFitParticle*, TVector3>& locDecayingParticleTrackPointGuesses);
		double Calc_TimeGuess(const deque<const DKinFitParticle*>& locFinalKinFitParticles, DVector3 locVertexGuess, bool locUseRFTimeFlag, double locRFTime);
		void Remove_BadVertexConstraints(deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_Vertices, deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_Vertices) const;

		void Build_KinFitResults(const DParticleCombo* locParticleCombo, const map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locInitDecayingKinFitParticles);

		const DAnalysisUtilities* dAnalysisUtilities;
		DKinFitter_GlueX dKinFitter;

		unsigned int dDebugLevel;
		unsigned int dKinFitDebugLevel;
		bool dLinkVerticesFlag;
};

#endif // _DKinFitResults_factory_

