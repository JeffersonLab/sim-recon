#ifndef _DKinFitResults_factory_
#define _DKinFitResults_factory_

#include <deque>
#include <vector>
#include <map>
#include <set>

#include "TString.h"

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

		//GET CURRENT POOL SIZES
		size_t Get_KinFitParticlePoolSize(void) const{return dKinFitter.Get_KinFitParticlePoolSize();};
		size_t Get_KinFitConstraintVertexPoolSize(void) const{return dKinFitter.Get_KinFitConstraintVertexPoolSize();};
		size_t Get_KinFitConstraintSpacetimePoolSize(void) const{return dKinFitter.Get_KinFitConstraintSpacetimePoolSize();};
		size_t Get_KinFitConstraintP4PoolSize(void) const{return dKinFitter.Get_KinFitConstraintP4PoolSize();};
		size_t Get_MatrixDSymPoolSize(void) const{return dKinFitter.Get_MatrixDSymPoolSize();};
		size_t Get_LargeMatrixDSymPoolSize(void) const{return dKinFitter.Get_LargeMatrixDSymPoolSize();};

	private:

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		bool Create_KinFitConstraints(const DParticleCombo* locParticleCombo, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locDecayingKinFitParticles, deque<DKinFitConstraint*>& locOriginalConstraints, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints);
		void Setup_P4Constraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_P4, deque<const DKinFitParticle*>& locFinalKinFitParticles_P4, deque<size_t>& locIncludedStepIndices, bool& locConstrainMassFlag);
		void Setup_VertexConstraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_Vertex, deque<const DKinFitParticle*>& locFinalKinFitParticles_Vertex, deque<size_t>& locIncludedStepIndices);

		bool Handle_IfKinFitResultsWillBeIdentical(const DParticleCombo* locParticleCombo, deque<DKinFitConstraint*> locConstraints_ToCheck, const DEventRFBunch* locRFBunch_ToCheck, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles_ToCheck);
		bool Check_IfKinFitResultsWillBeIdentical(map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles_ToCheck, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles_CheckAgainst);
		bool Check_IfKinFitResultsWillBeIdentical(deque<DKinFitConstraint*> locConstraints_ToCheck, deque<const DKinFitConstraint*> locConstraints_CheckAgainst, const DEventRFBunch* locRFBunch_ToCheck, const DEventRFBunch* locRFBunch_CheckAgainst);
		bool Check_IfKinFitResultsWillBeIdentical(DKinFitConstraint_P4* locConstraint_ToCheck, const DKinFitConstraint_P4* locConstraint_CheckAgainst);
		bool Check_IfKinFitResultsWillBeIdentical(DKinFitConstraint_VertexBase* locConstraint_ToCheck, const DKinFitConstraint_VertexBase* locConstraint_CheckAgainst);
		bool Check_IfKinFitResultsWillBeIdentical(DKinFitConstraint_Spacetime* locConstraint_ToCheck, const DKinFitConstraint_Spacetime* locConstraint_CheckAgainst, const DEventRFBunch* locRFBunch_ToCheck, const DEventRFBunch* locRFBunch_CheckAgainst);
		bool Check_IfKinFitResultsWillBeIdentical(deque<const DKinFitParticle*> locParticles_ToCheck, deque<const DKinFitParticle*> locParticles_CheckAgainst);

		bool Setup_KinFit(DKinFitType locKinFitType, const deque<DKinFitConstraint*>& locOriginalConstraints, const DEventRFBunch* locEventRFBunch, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints);
		double Calc_TimeGuess(const DKinFitConstraint_Spacetime* locConstraint, DVector3 locVertexGuess, double locRFTime);
		void Build_KinFitResults(const DParticleCombo* locParticleCombo, const map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locInitDecayingKinFitParticles, deque<DKinFitConstraint*>& locOriginalConstraints);

		const DAnalysisUtilities* dAnalysisUtilities;
		DKinFitter_GlueX dKinFitter;

		unsigned int dDebugLevel;
		unsigned int dKinFitDebugLevel;
		bool dLinkVerticesFlag;

		double dTargetZCenter;

		class DPreviousFitInfo
		{
			public:
				DPreviousFitInfo(const DEventRFBunch* locEventRFBunch, deque<const DKinFitConstraint*>& locOriginalConstraints, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locDecayingParticles) : 
				dEventRFBunch(locEventRFBunch), dOriginalConstraints(locOriginalConstraints), dDecayingParticles(locDecayingParticles) {}

				const DEventRFBunch* dEventRFBunch;
				deque<const DKinFitConstraint*> dOriginalConstraints; //post skim after sort
				map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > dDecayingParticles;
			private:
				DPreviousFitInfo(void);
		};

		deque<DPreviousFitInfo> dPreviouslyFailedFits;
};

#endif // _DKinFitResults_factory_

