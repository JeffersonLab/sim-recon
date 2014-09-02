// $Id: ANALYSIS_init.cc 2433 2007-04-07 14:57:32Z pmatt $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DReaction_factory_Thrown.h"
#include "DParticleCombo_factory_Thrown.h"

#include "DAnalysisUtilities_factory.h"
#include "DParticleComboBlueprint_factory.h"

#include "DTrackTimeBased_factory_Combo.h"
#include "DDetectorMatches_factory_Combo.h"
#include "DEventRFBunch_factory_Combo.h"
#include "DChargedTrackHypothesis_factory_Combo.h"
#include "DNeutralParticleHypothesis_factory_Combo.h"

#include "DParticleCombo_factory_PreKinFit.h"
#include "DKinFitResults_factory.h"

#include "DBeamPhoton_factory_KinFit.h"
#include "DChargedTrackHypothesis_factory_KinFit.h"
#include "DNeutralParticleHypothesis_factory_KinFit.h"

#include "DParticleCombo_factory.h"
#include "DMCThrownMatching_factory.h"
#include "DAnalysisResults_factory_PreKinFit.h"
#include "DAnalysisResults_factory.h"
#include "DEventWriterROOT_factory.h"

#include "DHistogramActions.h"
#include "DCutActions.h"

jerror_t ANALYSIS_init(JEventLoop *loop)
{
	/// Create and register ANALYSIS data factories
	loop->AddFactory(new DReaction_factory_Thrown);
	loop->AddFactory(new DParticleCombo_factory_Thrown);

	loop->AddFactory(new DAnalysisUtilities_factory);
	loop->AddFactory(new DParticleComboBlueprint_factory);

	loop->AddFactory(new DEventRFBunch_factory_Combo);
	loop->AddFactory(new DTrackTimeBased_factory_Combo);
	loop->AddFactory(new DDetectorMatches_factory_Combo);
	loop->AddFactory(new DChargedTrackHypothesis_factory_Combo);
	loop->AddFactory(new DNeutralParticleHypothesis_factory_Combo);

	loop->AddFactory(new DParticleCombo_factory_PreKinFit);
	loop->AddFactory(new DKinFitResults_factory);

	loop->AddFactory(new DBeamPhoton_factory_KinFit);
	loop->AddFactory(new DChargedTrackHypothesis_factory_KinFit);
	loop->AddFactory(new DNeutralParticleHypothesis_factory_KinFit);

	loop->AddFactory(new DParticleCombo_factory);
	loop->AddFactory(new DMCThrownMatching_factory);
	loop->AddFactory(new DAnalysisResults_factory_PreKinFit);
	loop->AddFactory(new DAnalysisResults_factory);
	loop->AddFactory(new DEventWriterROOT_factory);

	//For some reason, have difficulty linking these classes without using them somewhere within the library
	DHistogramAction_ThrownParticleKinematics();
	DHistogramAction_DetectedParticleKinematics();
	DHistogramAction_ReconnedThrownKinematics();
	DHistogramAction_GenReconTrackComparison();
	DHistogramAction_TrackMultiplicity();
	DHistogramAction_TOFHitStudy();
	DHistogramAction_NumReconstructedObjects();

	DHistogramAction_PID(NULL);
	DHistogramAction_TrackVertexComparison(NULL);
	DHistogramAction_ParticleComboKinematics(NULL, false);
	DHistogramAction_TruePID(NULL);
	DHistogramAction_InvariantMass(NULL, Unknown, false, 0, 0.0, 0.0);
	DHistogramAction_MissingMass(NULL, false, 0, 0.0, 0.0);
	DHistogramAction_MissingMassSquared(NULL, false, 0, 0.0, 0.0);
	DHistogramAction_KinFitResults(NULL, 0.0);
	DHistogramAction_NumParticleCombos(NULL);
	DHistogramAction_ParticleComboGenReconComparison(NULL, false);

	DCutAction_ThrownTopology(NULL, true);
	DCutAction_PIDFOM(NULL, Unknown, Unknown, 0.0);
	DCutAction_AllTracksHaveDetectorMatch(NULL);
	DCutAction_CombinedPIDFOM(NULL, 0.0);
	DCutAction_CombinedTrackingFOM(NULL, 0.0);
	DCutAction_MissingMass(NULL, false, 0.0, 0.0);
	DCutAction_MissingMassSquared(NULL, false, 0.0, 0.0);
	DCutAction_InvariantMass(NULL, Unknown, false, 0.0, 0.0);
	DCutAction_AllVertexZ(NULL, 0.0, 0.0);
	DCutAction_MaxTrackDOCA(NULL, Unknown, 0.0);
	DCutAction_KinFitFOM(NULL, 0.0);
	DCutAction_TruePID(NULL, Unknown, Unknown, 0.0);
	DCutAction_AllTruePID(NULL, 0.0);
	DCutAction_GoodEventRFBunch(NULL, false);
	DCutAction_TransverseMomentum(NULL, 0.0);

	return NOERROR;
}


