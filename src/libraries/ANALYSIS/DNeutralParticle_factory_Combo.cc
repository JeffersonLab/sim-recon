#include "DNeutralParticle_factory_Combo.h"

//------------------
// init
//------------------
jerror_t DNeutralParticle_factory_Combo::init(void)
{
	//Get preselect tag
	dShowerSelectionTag = "PreSelect";
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticle_factory_Combo::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses); //make sure that brun() is called for the default factory!!!
	dNeutralParticleHypothesisFactory = static_cast<DNeutralParticleHypothesis_factory*>(locEventLoop->GetFactory("DNeutralParticleHypothesis"));

	//Get Needed PIDs
	auto locReactions = DAnalysis::Get_Reactions(locEventLoop);
	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		auto locNeutralPIDs = locReactions[loc_i]->Get_FinalPIDs(-1, false, false, d_Neutral, false);
		for(auto locPID : locNeutralPIDs)
		{
			if(locPID != Gamma) //already created by default!
				dNeutralPIDs.insert(locPID);
		}
	}

	//Setting this flag makes it so that JANA does not delete the objects in _data.  This factory will manage this memory.
	if(dNeutralPIDs.empty())
		SetFactoryFlag(NOT_OBJECT_OWNER);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticle_factory_Combo::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, dShowerSelectionTag.c_str());

	//Nothing to do! Pass them through
	if(dNeutralPIDs.empty())
	{
		_data.clear();
		for(auto& locNeutralParticle : locNeutralParticles)
			_data.push_back(const_cast<DNeutralParticle*>(locNeutralParticle));
		return NOERROR;
	}

	const DEventRFBunch* locEventRFBunch = nullptr;
	locEventLoop->GetSingle(locEventRFBunch);

	const DVertex* locVertex = nullptr;
	locEventLoop->GetSingle(locVertex);

	dNeutralParticleHypothesisFactory->Recycle_Hypotheses(dCreatedHypotheses);

	//create and add new hypos
	for(auto& locNeutralParticle : locNeutralParticles)
	{
		auto locNewNeutralParticle = new DNeutralParticle(*locNeutralParticle);
		for(auto& locPID : dNeutralPIDs)
		{
			//create new DNeutralParticleHypothesis object
			auto locNewHypothesis = dNeutralParticleHypothesisFactory->Create_DNeutralParticleHypothesis(locNeutralParticle->dNeutralShower, locPID, locEventRFBunch, locVertex->dSpacetimeVertex, &locVertex->dCovarianceMatrix);
			dCreatedHypotheses.push_back(locNewHypothesis);
			locNewNeutralParticle->dNeutralParticleHypotheses.push_back(locNewHypothesis);
		}
		_data.push_back(locNewNeutralParticle);
	}

	return NOERROR;
}
