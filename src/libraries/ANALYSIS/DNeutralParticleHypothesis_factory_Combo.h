// $Id$
//
//    File: DNeutralParticleHypothesis_factory_Combo.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralParticleHypothesis_factory_Combo_
#define _DNeutralParticleHypothesis_factory_Combo_

#include <limits>

#include <JANA/JFactory.h>
#include <PID/DNeutralParticleHypothesis.h>
#include "PID/DNeutralParticleHypothesis_factory.h"
#include <PID/DNeutralShower.h>
#include <PID/DNeutralParticle.h>
#include <PID/DEventRFBunch.h>

#include "ANALYSIS/DParticleComboBlueprint.h"

class DNeutralParticleHypothesis_factory_Combo:public jana::JFactory<DNeutralParticleHypothesis>
{
	public:
		DNeutralParticleHypothesis_factory_Combo(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DNeutralParticleHypothesis_factory_Combo(){};
		const char* Tag(void){return "Combo";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void Build_RFPIDMap(map<Particle_t, vector<const DReaction*> >& locPIDMap, map<const DReaction*, set<const DEventRFBunch*> >& locRFBunchReactionMap, map<Particle_t, set<const DEventRFBunch*> >& locRFBunchPIDMap);

		DNeutralParticleHypothesis_factory* dNeutralParticleHypothesisFactory;

		vector<const DReaction*> dReactions;
		map<Particle_t, vector<const DReaction*> > dNeutralPIDs; //vector: reactions for which they are needed
		string dShowerSelectionTag;
};

#endif // _DNeutralParticleHypothesis_factory_Combo_

