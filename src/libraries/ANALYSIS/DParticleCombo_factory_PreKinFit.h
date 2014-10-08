// $Id$
//
//    File: DParticleCombo_factory_PreKinFit.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleCombo_factory_PreKinFit_
#define _DParticleCombo_factory_PreKinFit_

#include <deque>
#include <vector>
#include <map>

#include "JANA/JFactory.h"
#include "JANA/JObject.h"

#include "particleType.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DKinematicData.h"
#include "PID/DChargedTrack.h"
#include "PID/DBeamPhoton.h"
#include "PID/DNeutralShower.h"
#include "PID/DEventRFBunch.h"
#include "PID/DVertex.h"
#include <HDGEOMETRY/DGeometry.h>
#include <DANA/DApplication.h>

#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DCutActions.h"
#include "ANALYSIS/DMCThrownMatching.h"
#include "ANALYSIS/DParticleComboStep.h"
#include "ANALYSIS/DParticleComboBlueprint.h"
#include "ANALYSIS/DAnalysisUtilities.h"

class DParticleCombo_factory_PreKinFit : public jana::JFactory<DParticleCombo>
{
	public:
		DParticleCombo_factory_PreKinFit(){};
		~DParticleCombo_factory_PreKinFit(){};
		const char* Tag(void){return "PreKinFit";}

		void Reset_Pools(void);

		size_t Get_ParticleComboStepPoolSize(void) const{return dParticleComboStepPool_All.size();};
		size_t Get_KinematicDataPoolSize(void) const{return dKinematicDataPool_All.size();};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DKinematicData* Get_DetectedParticle(const DReaction* locReaction, const DEventRFBunch* locEventRFBunch, const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locParticleIndex);
		DKinematicData* Create_Target(Particle_t locPID);

		bool Cut_PIDFOM(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const;
		bool Cut_PIDFOM(const DReaction* locReaction, const DNeutralParticleHypothesis* locNeutralParticleHypothesis) const;

		void Build_BeamPhotonCombos(DParticleCombo* locParticleCombo, const DParticleComboBlueprint* locParticleComboBlueprint, const DEventRFBunch* locEventRFBunch, const set<const DBeamPhoton*>& locInputCandidatePhotons, vector<DParticleCombo*>& locBuiltParticleCombos);

		DParticleComboStep* Clone_ParticleComboStep(const DParticleComboStep* locParticleComboStep);
		void Reset_KinematicData(DKinematicData* locKinematicData);
		DParticleComboStep* Get_ParticleComboStepResource(void);
		DParticleCombo* Get_ParticleComboResource(void);
		DKinematicData* Get_KinematicDataResource(void);
		DBeamPhoton* Get_BeamPhotonResource(void);

		void Recycle_Data(const DParticleCombo* locParticleCombo, const DParticleComboBlueprint* locParticleComboBlueprint, bool locAllButFirstStepFlag);
		void Recycle_Data_BeamStep(const DParticleCombo* locParticleCombo, const DParticleComboBlueprint* locParticleComboBlueprint, const DEventRFBunch* locEventRFBunch);

		deque<DParticleComboStep*> dParticleComboStepPool_All;
		deque<DParticleComboStep*> dParticleComboStepPool_Available;

		deque<DKinematicData*> dKinematicDataPool_All;
		deque<DKinematicData*> dKinematicDataPool_Available;

		map<const DParticleComboBlueprintStep*, const DParticleComboStep*> dComboBlueprintStepMap;
		map<pair<const DParticleComboBlueprintStep*, const DEventRFBunch*>, set<const DParticleComboStep*> > dComboBlueprintBeamStepMap;

		size_t MAX_DParticleComboStepPoolSize;
		size_t MAX_DKinematicDataPoolSize;

		// PRE-DPARTICLECOMBO CUT VALUES
			//bool = true/false for cut enabled/disabled, double = cut value
			//Command-line values will override these values
		pair<bool, double> dMinChargedPIDFOM; //the minimum PID FOM for a particle used for this DReaction
		pair<bool, double> dMinPhotonPIDFOM; //the minimum PID FOM for a neutral particle used for this DReaction
		pair<bool, double> dMaxPhotonRFDeltaT; //the maximum photon-rf time difference: used for photon selection

		int dDebugLevel;
		double dTargetCenterZ;
		double dMinThrownMatchFOM;
		const DAnalysisUtilities* dAnalysisUtilities;

		vector<const DReaction*> dReactions;
		map<const DReaction*, bool> dMCReactionExactMatchFlags;
		map<const DReaction*, DCutAction_TrueCombo*> dTrueComboCuts;
		map<const DReaction*, size_t> dNumGoodPreComboSelectionActions;

		map<const DReaction*, map<const DEventRFBunch*, map<Particle_t, map<const DChargedTrack*, const DChargedTrackHypothesis*> > > > dValidChargedHypotheses;
		map<const DReaction*, map<const DEventRFBunch*, map<Particle_t, map<const DNeutralShower*, const DNeutralParticleHypothesis*> > > > dValidNeutralHypotheses;

		map<const DReaction*, map<Particle_t, TH1I*> > dHistMap_PIDFOM_All;
		map<const DReaction*, map<Particle_t, TH1I*> > dHistMap_PIDFOM_True;

		map<const DReaction*, TH1I*> dHistMap_PhotonRFDeltaT_All;
		map<const DReaction*, TH1I*> dHistMap_PhotonRFDeltaT_True; //for the true combo (MC event) (includes true beam photon): delta-t between selected photon & RF time

		map<const DReaction*, TH1D*> dHistMap_NumSurvivingBeamParticles;

		map<const DReaction*, TH2D*> dHistMap_NumBlueprintsSurvivedCut;
		map<const DReaction*, TH1D*> dHistMap_NumBlueprintsSurvivedCut1D;

		map<const DReaction*, TH1D*> dHistMap_NumEventsSurvivedCut_All;
		map<const DReaction*, TH1D*> dHistMap_NumEventsSurvivedCut_True;
};

#endif // _DParticleCombo_factory_PreKinFit_

