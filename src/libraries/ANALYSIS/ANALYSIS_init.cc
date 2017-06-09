// $Id: ANALYSIS_init.cc 2433 2007-04-07 14:57:32Z pmatt $

#include <JANA/JEventLoop.h>
using namespace jana;

//OK
#include "DReaction_factory_Thrown.h"
#include "DParticleCombo_factory_Thrown.h"

#include "DAnalysisUtilities_factory.h"
#include "DMCThrownMatching_factory.h"

#include "DAnalysisResults_factory.h"
#include "DEventWriterROOT_factory.h"

#include "DTrackTimeBased_factory_Combo.h"
#include "DDetectorMatches_factory_Combo.h"
#include "DChargedTrack_factory_Combo.h"
#include "DNeutralParticle_factory_Combo.h"

#include "DReactionVertexInfo_factory.h"

#include "DHistogramActions.h"
#include "DCutActions.h"

jerror_t ANALYSIS_init(JEventLoop *loop)
{
	/// Create and register ANALYSIS data factories
	loop->AddFactory(new DReaction_factory_Thrown);
	loop->AddFactory(new DParticleCombo_factory_Thrown);

	loop->AddFactory(new DAnalysisUtilities_factory);
	loop->AddFactory(new DMCThrownMatching_factory);

	loop->AddFactory(new DAnalysisResults_factory);
	loop->AddFactory(new DEventWriterROOT_factory);

	loop->AddFactory(new DTrackTimeBased_factory_Combo);
	loop->AddFactory(new DDetectorMatches_factory_Combo);
	loop->AddFactory(new DChargedTrack_factory_Combo);
	loop->AddFactory(new DNeutralParticle_factory_Combo);

	loop->AddFactory(new DReactionVertexInfo_factory);

	//For some reason, have difficulty linking these classes without using them somewhere within the library
	DHistogramAction_ThrownParticleKinematics();
	DHistogramAction_DetectedParticleKinematics();
	DHistogramAction_ReconnedThrownKinematics();
	DHistogramAction_GenReconTrackComparison();
	DHistogramAction_TrackMultiplicity();
	DHistogramAction_TOFHitStudy();
	DHistogramAction_NumReconstructedObjects();
	DHistogramAction_DetectorMatchParams();
	DHistogramAction_Neutrals();
	DHistogramAction_DetectorPID();
	DHistogramAction_DetectorMatching();
	DHistogramAction_Reconstruction();
	DHistogramAction_ObjectMemory();

	DHistogramAction_PID(NULL);
	DHistogramAction_TrackVertexComparison(NULL);
	DHistogramAction_ParticleComboKinematics(NULL, false);
	DHistogramAction_TruePID(NULL);
	DHistogramAction_InvariantMass(NULL, Unknown, false, 0, 0.0, 0.0);
	DHistogramAction_MissingMass(NULL, false, 0, 0.0, 0.0);
	DHistogramAction_MissingMassSquared(NULL, false, 0, 0.0, 0.0);
	DHistogramAction_KinFitResults(NULL, 0.0);
	DHistogramAction_ParticleComboGenReconComparison(NULL, false);
	DHistogramAction_MissingTransverseMomentum(NULL, false, 0, 0.0, 0.0);
	DHistogramAction_2DInvariantMass(NULL, 0, deque<Particle_t>(), deque<Particle_t>(), false, 0, 0.0, 0.0, 0, 0.0, 0.0);
	DHistogramAction_Dalitz(NULL, 0, deque<Particle_t>(), deque<Particle_t>(), false, 0, 0.0, 0.0, 0, 0.0, 0.0);

	DCutAction_MinTrackHits(NULL, 0);
	DCutAction_ThrownTopology(NULL, true);
	DCutAction_PIDFOM(NULL, Unknown, Unknown, 0.0);
	DCutAction_AllTracksHaveDetectorMatch(NULL);
	DCutAction_CombinedPIDFOM(NULL, 0.0);
	DCutAction_EachPIDFOM(NULL, 0.0);
	DCutAction_CombinedTrackingFOM(NULL, 0.0);
	DCutAction_MissingMass(NULL, false, 0.0, 0.0);
	DCutAction_MissingMassSquared(NULL, false, 0.0, 0.0);
	DCutAction_InvariantMass(NULL, Unknown, false, 0.0, 0.0);
	DCutAction_AllVertexZ(NULL, 0.0, 0.0);
	DCutAction_ProductionVertexZ(NULL, 0.0, 0.0);
	DCutAction_MaxTrackDOCA(NULL, Unknown, 0.0);
	DCutAction_KinFitFOM(NULL, 0.0);
	DCutAction_TruePID(NULL, Unknown, Unknown, 0.0);
	DCutAction_AllTruePID(NULL, 0.0);
	DCutAction_GoodEventRFBunch(NULL, false);
	DCutAction_TransverseMomentum(NULL, 0.0);
	DCutAction_TrueBeamParticle(NULL);
	DCutAction_TrueCombo(NULL, 0.0, false);
	DCutAction_BDTSignalCombo(NULL, 0.0, false, false);

	DCutAction_TrackHitPattern(NULL);
	DCutAction_ProtonPiPlusdEdx(NULL, 0.0);
	DCutAction_BeamEnergy(NULL, false, 0.0, 0.0);
	DCutAction_TrackFCALShowerEOverP(NULL, false, 0.0);
	DCutAction_NoPIDHit(NULL);
	DCutAction_PIDDeltaT(NULL, false, 0.0);
	DCutAction_PIDTimingBeta(NULL, 0.0, 0.0);
	DCutAction_OneVertexKinFit(NULL);

	return NOERROR;
}

