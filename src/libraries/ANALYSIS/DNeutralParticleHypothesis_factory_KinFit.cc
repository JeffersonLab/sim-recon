// $Id$
//
//    File: DNeutralParticleHypothesis_factory_KinFit.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DNeutralParticleHypothesis_factory_KinFit.h"

//------------------
// init
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	locEventLoop->GetSingle(dParticleID);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DNeutralParticleHypothesis_factory_KinFit::evnt()");
#endif

 	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	map<DKinFitParticle*, DNeutralParticleHypothesis*> locNewObjectMap;

	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		map<const DParticleCombo*, const DKinFitChain*> locParticleComboMap;
		locKinFitResultsVector[loc_i]->Get_ParticleComboMap(locParticleComboMap);
		set<DKinFitParticle*> locOutputKinFitParticles = locKinFitResultsVector[loc_i]->Get_OutputKinFitParticles();

		map<const DParticleCombo*, const DKinFitChain*>::iterator locComboIterator = locParticleComboMap.begin();
		for(; locComboIterator != locParticleComboMap.end(); ++locComboIterator)
		{
			const DParticleCombo* locParticleCombo = locComboIterator->first;
			for(size_t loc_j = 0; loc_j < locParticleCombo->Get_NumParticleComboSteps(); ++loc_j)
			{
				const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_j);
				for(size_t loc_k = 0; loc_k < locParticleComboStep->Get_NumFinalParticles(); ++loc_k)
				{
					if(!locParticleComboStep->Is_FinalParticleDetected(loc_k))
						continue;
					if(!locParticleComboStep->Is_FinalParticleNeutral(loc_k))
						continue;

					//might have used neutral shower OR neutral particle. try particle first
					const DNeutralShower* locNeutralShower = static_cast<const DNeutralShower*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_k));
					const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticleComboStep->Get_FinalParticle(loc_k));

					DKinFitParticle* locKinFitParticle = locKinFitResultsVector[loc_i]->Get_OutputKinFitParticle(locNeutralParticleHypothesis);
					if(locKinFitParticle == NULL)
						locKinFitParticle = locKinFitResultsVector[loc_i]->Get_OutputKinFitParticle(locNeutralShower);
					if(locKinFitParticle == NULL)
						continue; //should be impossible
					if(locOutputKinFitParticles.find(locKinFitParticle) == locOutputKinFitParticles.end())
						continue; //not used in fit

					map<DKinFitParticle*, DNeutralParticleHypothesis*>::iterator locNewHypoIterator = locNewObjectMap.find(locKinFitParticle);
					if(locNewHypoIterator != locNewObjectMap.end())
					{
						locNewHypoIterator->second->AddAssociatedObject(locParticleCombo);
						continue; //new particle already created for this kinfit particle
					}

					DNeutralParticleHypothesis* locNewNeutralParticleHypothesis = Build_NeutralParticleHypothesis(locNeutralParticleHypothesis, locKinFitParticle, locNeutralShower, locParticleCombo);
					locNewObjectMap[locKinFitParticle] = locNewNeutralParticleHypothesis;
					locNewNeutralParticleHypothesis->AddAssociatedObject(locParticleCombo);

					_data.push_back(locNewNeutralParticleHypothesis);
				}
			}
		}
	}

	return NOERROR;
}

DNeutralParticleHypothesis* DNeutralParticleHypothesis_factory_KinFit::Build_NeutralParticleHypothesis(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, DKinFitParticle* locKinFitParticle, const DNeutralShower* locNeutralShower, const DParticleCombo* locParticleCombo)
{
	DNeutralParticleHypothesis* locNewNeutralParticleHypothesis = new DNeutralParticleHypothesis(*locNeutralParticleHypothesis);
	locNewNeutralParticleHypothesis->AddAssociatedObject(locNeutralParticleHypothesis);
	locNewNeutralParticleHypothesis->AddAssociatedObject(locNeutralShower);

 	vector<const JObject*> locObjects;
	locNeutralParticleHypothesis->GetT(locObjects);
	for(size_t loc_i = 0; loc_i < locObjects.size(); ++loc_i)
	{
		if(dynamic_cast<const DNeutralParticleHypothesis*>(locObjects[loc_i]) != NULL)
			continue; //don't save: won't be able to keep track of which is which! (combo or default factory)
		locNewNeutralParticleHypothesis->AddAssociatedObject(locObjects[loc_i]);
	}

	locNewNeutralParticleHypothesis->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewNeutralParticleHypothesis->setPosition(DVector3(locKinFitParticle->Get_CommonVertex().X(),locKinFitParticle->Get_CommonVertex().Y(),locKinFitParticle->Get_CommonVertex().Z()));
	locNewNeutralParticleHypothesis->setTime(locKinFitParticle->Get_CommonTime());
	locNewNeutralParticleHypothesis->setErrorMatrix(*locKinFitParticle->Get_CovarianceMatrix());

	if(locKinFitParticle->Get_ShowerEnergy() > 0.0) //particle was used in the fit as a neutral shower
		locNewNeutralParticleHypothesis->setPathLength(locKinFitParticle->Get_PathLength(), locKinFitParticle->Get_PathLengthUncertainty());
	else
	{
		double locPathLength = locNewNeutralParticleHypothesis->pathLength() - locKinFitParticle->Get_PathLength();
		double locPathLengthUncertainty_Orig = locNewNeutralParticleHypothesis->pathLength_err();
		double locPathLengthUncertainty_KinFit = locKinFitParticle->Get_PathLengthUncertainty();
		double locPathLengthUncertainty = sqrt(locPathLengthUncertainty_Orig*locPathLengthUncertainty_Orig + locPathLengthUncertainty_KinFit*locPathLengthUncertainty_KinFit);
		locNewNeutralParticleHypothesis->setPathLength(locPathLength, locPathLengthUncertainty);
	}

	// Calculate DNeutralParticleHypothesis FOM
	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();
	double locPropagatedRFTime = dParticleID->Calc_PropagatedRFTime(locNewNeutralParticleHypothesis, locEventRFBunch);
	locNewNeutralParticleHypothesis->setT0(locPropagatedRFTime, sqrt(locEventRFBunch->dTimeVariance), locEventRFBunch->dTimeSource);

	// Calculate DNeutralParticleHypothesis FOM
	unsigned int locNDF = 0;
	double locChiSq = 0.0;
	double locFOM = -1.0; //undefined for non-photons
	if(locNewNeutralParticleHypothesis->PID() == Gamma)
	{
		double locTimePull = 0.0;
		//for this calc: if rf time part of timing constraint, don't use locKinFitParticle->Get_CommonTime() for chisq calc!!!
		locChiSq = dParticleID->Calc_TimingChiSq(locNewNeutralParticleHypothesis, locNDF, locTimePull);
		locFOM = TMath::Prob(locChiSq, locNDF);
	}
	locNewNeutralParticleHypothesis->dChiSq = locChiSq;
	locNewNeutralParticleHypothesis->dNDF = locNDF;
	locNewNeutralParticleHypothesis->dFOM = locFOM;

	return locNewNeutralParticleHypothesis;
}

//------------------
// erun
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::fini(void)
{
	return NOERROR;
}


